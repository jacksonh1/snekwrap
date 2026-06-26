"""
Thin functional wrapper around PottsMPNN energy prediction (analogous to
``run_potts_mpnn``).

``run_energy_prediction`` takes an input PDB and all configuration as plain function
arguments, runs PottsMPNN's ``energy_prediction`` (a deep-mutational-scan saturation
sweep by default), and returns an :class:`EnergyPredictionResult`. The result holds the
per-mutation score table as a DataFrame and can draw the per-position x per-amino-acid
heatmap on demand via :meth:`EnergyPredictionResult.plot_heatmap`.

The PottsMPNN science code is NOT vendored here; it lives in the repo/submodule
(``config.POTTSMPNN_REPO``) and is imported lazily inside the context managers below,
exactly like ``potts_mpnn.py``.
"""
from __future__ import annotations

import tempfile
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from omegaconf import OmegaConf

import snekwrap.config as config

# Reuse the sys.path / cwd helpers from the sequence-design wrapper rather than
# duplicating them. Importing potts_mpnn does not import this module, so there is no
# circular import.
from snekwrap.wrappers.pottsmpnn.potts_mpnn import _temporary_sys_path, _temporary_cwd


_SCORE_COLUMNS = ["pdb", "mutant", "wildtype", "ddG_pred", "ddG_expt"]

# The 20 standard amino acids, in the same order PottsMPNN's run_utils.plot_data uses for
# the heatmap rows. _heatmap_grid() reproduces that ordering so the interactive heatmap
# matches the matplotlib one row-for-row.
_HEATMAP_AMINO_ACIDS = [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
]


@dataclass
class EnergyPredictionResult:
    """PottsMPNN energy-prediction results for a single PDB.

    The raw mutation scores live in ``scores`` (one row per scored mutant, as produced by
    PottsMPNN). For analysis, :attr:`mutation_energies` exposes the same predictions as a tidy
    DataFrame per chain (``position, wildtype, mutant, ddG_pred``). Use :meth:`plot_heatmap` or
    :meth:`plot_heatmap_interactive` to render the per-position x per-amino-acid heatmap.
    """

    pdb_name: str
    scores: pd.DataFrame  # columns: pdb, mutant, wildtype, ddG_pred, ddG_expt
    chain_order: List[str]  # PDB chain order; maps colon-delimited sequences to chains
    wildtype_sequence: str  # colon-joined per-chain wild-type sequence
    ener_type: str  # "ddG" (relative, centered at 0) or "dG" (absolute)
    repo_root: Path  # PottsMPNN repo/submodule root, for lazily importing plot_data
    input_params: Dict
    stats: Optional[pd.DataFrame] = None  # per-PDB Pearson r, if experimental values given

    def plot_heatmap(
        self,
        save_path: str | Path | None = None,
        chain_ranges: Optional[Dict[str, List[int]]] = None,
        figsize: tuple = (20, 5),
        title: str = "PottsMPNN Predictions",
        clabel: str = r"Predicted $\Delta\Delta$G (a.u.)",
        only_mutated_positions: bool = False,
    ) -> "matplotlib.figure.Figure":
        """Draw the mutation-energy heatmap (positions x 20 amino acids).

        Calls PottsMPNN's ``plot_data`` on the stored score table and returns the matplotlib
        ``Figure`` so it can be customized further (titles, annotations, resaving, ...). The
        heatmap axes are ``fig.axes[0]``.

        In a notebook, returning the figure as the last expression in a cell displays it once;
        assign it (``fig = result.plot_heatmap()``) to suppress display, modify, then redisplay.

        Parameters
        ----------
        save_path:
            If given, the heatmap is also saved here (PNG).
        chain_ranges:
            Optional ``{chain: [start, stop]}`` (inclusive, 1-indexed) to zoom the heatmap
            to specific position ranges.
        only_mutated_positions:
            If True, only show columns (positions) that were actually scanned. Useful when a
            subset of a large complex is scanned (e.g. ``scan_chains``), to avoid hundreds
            of empty columns.

        Returns
        -------
        matplotlib.figure.Figure
        """
        import matplotlib.pyplot as plt

        with _temporary_sys_path(self.repo_root):
            from run_utils import plot_data

            # plot_data builds its own figure and shows/closes it internally (it does not
            # return one). With verbose=True it routes to plt.show() rather than plt.close(),
            # so we no-op plt.show for the duration and then grab the surviving figure.
            _orig_show = plt.show
            plt.show = lambda *args, **kwargs: None
            try:
                plot_data(
                    self.scores,
                    only_mutated_positions=only_mutated_positions,
                    title=title,
                    clabel=clabel,
                    save_path=str(save_path) if save_path is not None else None,
                    figsize=figsize,
                    ener_type=self.ener_type,
                    chain_ranges=chain_ranges,
                    chain_order=self.chain_order,
                    verbose=True,
                )
                fig = plt.gcf()
            finally:
                plt.show = _orig_show

        # Detach from pyplot so the inline backend does not auto-render it a second time; the
        # Figure object stays fully usable (edit axes, savefig, display).
        plt.close(fig)
        return fig

    def _parse_scores(self) -> tuple[Dict[str, Dict[int, dict]], Dict[str, List[str]]]:
        """Parse ``self.scores`` into per-chain, per-position mutations.

        Shared source of truth behind both :attr:`mutation_energies` (the tidy per-chain tables)
        and :meth:`_heatmap_grid` (the wide grid the heatmaps render).

        Returns ``(parsed, chain_sequences)`` where
        ``parsed[chain][pos] = {"wt": wt_aa, "muts": {mut_aa: energy}}`` (``pos`` is 1-indexed)
        and ``chain_sequences[chain]`` is that chain's wild-type residue list.
        """
        # parsed[chain][pos] = {"wt": wt_aa, "muts": {mut_aa: energy}}
        parsed: Dict[str, Dict[int, dict]] = {}
        chain_sequences: Dict[str, List[str]] = {}

        # NOTE: the skip semantics below (silently ignoring rows whose chain counts/lengths
        # mismatch or that aren't single mutations) intentionally mirror run_utils.plot_data so
        # the grids/tables built from this are identical to the matplotlib heatmap. This is a
        # deliberate, documented exception to the project's fail-loudly default, for visual parity.
        for row in self.scores.itertuples(index=False):
            wt_chains = row.wildtype.split(":")
            mut_chains = row.mutant.split(":")
            if len(wt_chains) != len(mut_chains):
                continue

            chain_names = self.chain_order[: len(wt_chains)]
            mutations = []  # (chain, 1-indexed pos, wt_aa, mut_aa)
            for c_name, w_chain, m_chain in zip(chain_names, wt_chains, mut_chains, strict=True):
                if len(w_chain) != len(m_chain):
                    break
                chain_sequences.setdefault(c_name, list(w_chain))
                for i, (w, m) in enumerate(zip(w_chain, m_chain, strict=True)):
                    if w != m:
                        mutations.append((c_name, i + 1, w, m))

            # Only single-substitution rows map onto the grid.
            if len(mutations) != 1:
                continue
            c_name, pos, wt_aa, mut_aa = mutations[0]
            parsed.setdefault(c_name, {}).setdefault(pos, {"wt": wt_aa, "muts": {}})
            parsed[c_name][pos]["muts"][mut_aa] = row.ddG_pred

        return parsed, chain_sequences

    @cached_property
    def mutation_energies(self) -> Dict[str, pd.DataFrame]:
        """Per-chain tidy table of the mutation-energy predictions.

        Lazily computed on first access and cached on the instance. Returns
        ``{chain: DataFrame}`` (chains in ``self.chain_order`` order, only those present in the
        data). Each DataFrame has one row per scored substitution with columns
        ``["position", "wildtype", "mutant", "ddG_pred"]``, sorted by position then mutant
        residue. Wild-type "self" rows are not included. Slice it as plain data, e.g.
        ``result.mutation_energies["A"].query("position == 42")`` or
        ``df[df.position.between(10, 20)]``.
        """
        aa_order = {aa: i for i, aa in enumerate(_HEATMAP_AMINO_ACIDS)}
        parsed, _ = self._parse_scores()
        # Only chains that actually have scored mutations get a table. (Every chain shows up in
        # `chain_sequences` because it appears in each row's unchanged sequence, but a chain that
        # was never scanned -- e.g. omitted via `scan_chains` -- has no mutations and would
        # otherwise yield an empty, noisy table.)
        scanned_chains = [c for c in self.chain_order if c in parsed]

        tables: Dict[str, pd.DataFrame] = {}
        for c_name in scanned_chains:
            rows = []
            for pos in sorted(parsed[c_name]):
                entry = parsed[c_name][pos]
                wt_aa = entry["wt"]
                for mut_aa, energy in sorted(entry["muts"].items(), key=lambda kv: aa_order[kv[0]]):
                    rows.append((pos, wt_aa, mut_aa, energy))
            tables[c_name] = pd.DataFrame(
                rows, columns=["position", "wildtype", "mutant", "ddG_pred"]
            )
        return tables

    def _heatmap_grid(
        self,
        chain_ranges: Optional[Dict[str, List[int]]] = None,
        only_mutated_positions: bool = False,
    ) -> pd.DataFrame:
        """Build the [20 amino acids x positions] grid the heatmaps render.

        Reproduces the parsing and column selection of PottsMPNN's ``run_utils.plot_data`` so the
        interactive heatmap matches the matplotlib one row-for-row. Index is the 20 amino acids
        (mutant residue); columns are a MultiIndex ``(chain, position, wildtype)``. Values are the
        predicted energy; NaN where a substitution was not scored. For ``ener_type == "ddG"`` the
        wild-type residue's cell in each column is 0.0.

        ``chain_ranges`` is ``{chain: [start, stop]}`` (inclusive, 1-indexed) restricting columns;
        ``start=0`` -> 1, ``stop=-1`` -> chain length; chains absent from the dict are dropped
        (matches ``plot_data``). ``only_mutated_positions`` keeps only scanned positions.
        """
        aa_to_idx = {aa: i for i, aa in enumerate(_HEATMAP_AMINO_ACIDS)}
        parsed, chain_sequences = self._parse_scores()

        # Order chains by chain_order, keeping only those that appeared in the data.
        active_chains = [c for c in self.chain_order if c in chain_sequences]

        # Build the ordered list of columns: (chain, pos, wt_aa).
        columns: List[tuple] = []
        for c_name in active_chains:
            full_seq = chain_sequences[c_name]
            start_r, stop_r = 1, len(full_seq)
            if chain_ranges and c_name in chain_ranges:
                start_r, stop_r = chain_ranges[c_name]
                if start_r == 0:
                    start_r = 1
                if stop_r == -1:
                    stop_r = len(full_seq)
            elif chain_ranges:
                # chain_ranges given but this chain not listed -> drop it (matches plot_data).
                continue

            if only_mutated_positions:
                positions = [p for p in sorted(parsed.get(c_name, {})) if start_r <= p <= stop_r]
            else:
                actual_start = max(1, start_r)
                actual_stop = min(len(full_seq), stop_r)
                positions = range(actual_start, actual_stop + 1) if actual_start <= actual_stop else []

            for pos in positions:
                columns.append((c_name, pos, full_seq[pos - 1]))

        grid = np.full((len(_HEATMAP_AMINO_ACIDS), len(columns)), np.nan)
        for col_idx, (c_name, pos, wt_aa) in enumerate(columns):
            if self.ener_type == "ddG" and wt_aa in aa_to_idx:
                grid[aa_to_idx[wt_aa], col_idx] = 0.0
            if c_name in parsed and pos in parsed[c_name]:
                for mut_aa, energy in parsed[c_name][pos]["muts"].items():
                    if mut_aa in aa_to_idx:
                        grid[aa_to_idx[mut_aa], col_idx] = energy

        col_index = pd.MultiIndex.from_tuples(columns, names=["chain", "position", "wildtype"])
        return pd.DataFrame(grid, index=list(_HEATMAP_AMINO_ACIDS), columns=col_index)

    def plot_heatmap_interactive(
        self,
        save_path: str | Path | None = None,
        chain_ranges: Optional[Dict[str, List[int]]] = None,
        title: str = "PottsMPNN Predictions",
        clabel: str = "Predicted ΔΔG (a.u.)",
        only_mutated_positions: bool = False,
        width: int | None = None,
        height: int = 500,
    ) -> "plotly.graph_objects.Figure":
        """Draw an interactive (plotly) mutation-energy heatmap and return the ``Figure``.

        Same data and layout as :meth:`plot_heatmap`, but interactive: hover a cell to read
        the exact predicted energy for that position/substitution. Returns a plotly
        ``graph_objects.Figure`` so you can ``.show()`` it, restyle it, or re-save it.

        Parameters
        ----------
        save_path:
            If given, the figure is written here. ``.html`` writes a standalone interactive
            file; any other extension uses ``fig.write_image`` (which requires ``kaleido``).
        chain_ranges:
            Optional ``{chain: [start, stop]}`` (inclusive, 1-indexed) to zoom to position
            ranges. ``start=0`` is treated as 1 and ``stop=-1`` as the chain length; chains absent
            from the dict are dropped.
        only_mutated_positions:
            If True, show only scanned positions (useful when scanning a subset of a complex).
        width, height:
            Figure size in pixels. ``width`` defaults to a value that scales with the number of
            columns so labels stay legible.

        Returns
        -------
        plotly.graph_objects.Figure
        """
        import plotly.graph_objects as go

        matrix = self._heatmap_grid(
            chain_ranges=chain_ranges, only_mutated_positions=only_mutated_positions
        )
        amino_acids = list(matrix.index)
        chains = matrix.columns.get_level_values("chain").tolist()
        positions = matrix.columns.get_level_values("position").tolist()
        wildtypes = matrix.columns.get_level_values("wildtype").tolist()
        n_cols = matrix.shape[1]
        z = matrix.values

        tick_labels = [f"{wt}{pos}" for wt, pos in zip(wildtypes, positions, strict=True)]

        # Per-cell hover metadata: (chain, position, wildtype, mutant_aa). customdata broadcasts
        # the per-column info across all 20 rows; the mutant residue is the row label.
        customdata = np.empty((len(amino_acids), n_cols, 4), dtype=object)
        for col_idx in range(n_cols):
            for row_idx, mut_aa in enumerate(amino_acids):
                customdata[row_idx, col_idx] = (
                    chains[col_idx], positions[col_idx], wildtypes[col_idx], mut_aa,
                )

        hovertemplate = (
            "Chain %{customdata[0]} · "
            "%{customdata[2]}%{customdata[1]}%{customdata[3]}<br>"
            f"{clabel}: " + "%{z:.2f}<extra></extra>"
        )

        # Diverging blue -> light gray -> red, matching run_utils.plot_data.
        colorscale = [[0.0, "blue"], [0.5, "rgb(230,230,230)"], [1.0, "red"]]
        zmid = 0.0 if self.ener_type == "ddG" else float(np.nanmean(z))

        fig = go.Figure(
            go.Heatmap(
                z=z,
                x=list(range(n_cols)),
                y=amino_acids,
                customdata=customdata,
                hovertemplate=hovertemplate,
                colorscale=colorscale,
                zmid=zmid,
                colorbar={"title": clabel},
                xgap=0,
                ygap=0,
            )
        )

        # Chain boundaries + labels (mirrors plot_data's per-chain borders).
        boundaries = [i for i in range(1, n_cols) if chains[i] != chains[i - 1]]
        for b in boundaries:
            fig.add_vline(x=b - 0.5, line={"color": "black", "width": 2})
        # Centered "Chain X" annotation per contiguous segment.
        seg_start = 0
        for end in boundaries + [n_cols]:
            c_name = chains[seg_start]
            fig.add_annotation(
                x=(seg_start + end - 1) / 2, y=1.02, yref="paper",
                text=f"Chain {c_name}", showarrow=False,
                font={"size": 13, "color": "black"},
            )
            seg_start = end

        fig.update_layout(
            title=title,
            width=width if width is not None else max(600, 18 * n_cols),
            height=height,
            plot_bgcolor="#E0E0E0",  # so NaN (un-scanned) cells read as light gray
            xaxis={
                "tickmode": "array",
                "tickvals": list(range(n_cols)),
                "ticktext": tick_labels,
                "tickangle": 90,
                "title": "Wildtype residue",
                "showgrid": False,
            },
            yaxis={"title": "Mutant residue", "autorange": "reversed", "showgrid": False},
        )

        if save_path is not None:
            save_path = str(save_path)
            if save_path.lower().endswith(".html"):
                fig.write_html(save_path)
            else:
                # write_image needs the optional `kaleido` package installed.
                fig.write_image(save_path)

        return fig

    def position_means(
        self,
        chains: Optional[List[str]] = None,
        agg: str = "mean",
    ) -> Dict[str, pd.Series]:
        """Summarize the per-position predicted energy for each chain.

        Aggregates the scored substitutions at every position down to a single number (the mean
        over the ~19 mutations by default), giving the per-residue signal behind
        :meth:`view_structure`. The reusable companion to :attr:`mutation_energies`.

        Parameters
        ----------
        chains:
            Restrict to these chains (default: all chains present in ``mutation_energies``).
        agg:
            How to collapse each position's substitutions: ``"mean"`` (default), ``"min"`` (the
            most-stabilizing substitution), ``"max"``, or any pandas-groupby aggregation
            string / callable.

        Returns
        -------
        dict[str, pandas.Series]
            ``{chain: Series}`` where each Series is indexed by position and holds the aggregated
            ``ddG_pred``.
        """
        me = self.mutation_energies
        chains = list(me) if chains is None else chains
        out: Dict[str, pd.Series] = {}
        for c in chains:
            if c not in me:
                raise ValueError(f"chain {c!r} has no scored mutations; available: {list(me)}.")
            out[c] = me[c].groupby("position")["ddG_pred"].agg(agg)
        return out

    def view_structure(
        self,
        chains: Optional[List[str]] = None,
        agg: str = "mean",
        cmap=None,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
        missing_color: str = "#bbbbbb",
        representation: str = "cartoon",
    ) -> "nglview.NGLWidget":
        """Color the input structure by per-position predicted energy and return an nglview widget.

        Each scanned residue is colored by its aggregated ``ddG_pred`` (mean over the saturation
        scan by default); residues that were never scored are ``missing_color`` (gray). The
        diverging scale is centered at 0 for ``ener_type == "ddG"`` (blue = stabilizing / negative,
        red = destabilizing / positive), else at the mean.

        The ``position`` index in the scores is PottsMPNN's 1-based sequence index; this method
        maps it to author residue numbers (``first_resnum + position - 1``) and **validates** that
        the structure's residue identity matches the recorded wild-type residue, raising if not.
        It therefore requires the run to have used ``skip_gaps=False`` (the default), so positions
        align to residue numbers.

        Parameters
        ----------
        chains:
            Restrict coloring to these chains (default: all scanned chains). Other chains/residues
            render in ``missing_color``.
        agg:
            Per-position aggregation passed to :meth:`position_means` (``"mean"``, ``"min"``, ...).
        cmap:
            A matplotlib colormap (name or object). Default: a blue->light-gray->red diverging map
            matching the heatmaps.
        vmin, vmax:
            Color limits. Default: symmetric about the center by the largest absolute deviation.
        missing_color:
            Color for residues with no predictions.
        representation:
            nglview representation for the colored structure (default ``"cartoon"``).

        Returns
        -------
        nglview.NGLWidget
        """
        import io

        import matplotlib.colors as mcolors
        import nglview as nv
        import biotite.structure as struc
        import biotite.structure.io.pdb as pdb_io
        from biotite.sequence import ProteinSequence

        if self.input_params["skip_gaps"]:
            raise ValueError(
                "view_structure requires the run to use skip_gaps=False so sequence positions map "
                "to author residue numbers; this result was produced with skip_gaps=True."
            )

        means = self.position_means(chains, agg)

        # Wild-type sequence per chain (skip_gaps=False -> '-' fills gaps, so length spans the full
        # residue-number range), aligned to chain_order.
        wt_by_chain = dict(zip(self.chain_order, self.wildtype_sequence.split(":"), strict=True))

        # Read the structure once; build per-chain residue-number info and a residue-identity map.
        pdb_path = self.input_params["pdb_path"]
        arr = pdb_io.PDBFile.read(pdb_path).get_structure(model=1)
        arr = arr[struc.filter_amino_acids(arr)]
        first_resnum: Dict[str, int] = {}
        resname_at: Dict[tuple, str] = {}  # (chain, res_id) -> 1-letter
        for chain_id, res_id, res_name, ins_code in zip(
            arr.chain_id, arr.res_id, arr.res_name, arr.ins_code, strict=True
        ):
            chain_id = str(chain_id)
            if chain_id not in means:
                continue
            if ins_code.strip():
                raise ValueError(
                    f"chain {chain_id!r} has insertion codes, which break the position->residue-"
                    "number mapping; view_structure does not support them."
                )
            key = (chain_id, int(res_id))
            if key not in resname_at:
                resname_at[key] = ProteinSequence.convert_letter_3to1(res_name)
            first_resnum[chain_id] = min(first_resnum.get(chain_id, int(res_id)), int(res_id))

        # Sanity check: PottsMPNN's gap-filled sequence length must span the residue-number range.
        for c in means:
            res_ids = [r for (ch, r) in resname_at if ch == c]
            span = max(res_ids) - first_resnum[c] + 1
            if span != len(wt_by_chain[c]):
                raise ValueError(
                    f"chain {c}: residue-number span ({span}) != wild-type sequence length "
                    f"({len(wt_by_chain[c])}); cannot map positions to residues reliably."
                )

        # Map each scored position to its residue number, validating that the structure's residue
        # identity matches the recorded wild type (catches a bad position->residue mapping).
        resno_value: Dict[tuple, float] = {}
        for c, series in means.items():
            wt_seq = wt_by_chain[c]
            for pos, value in series.items():
                resno = first_resnum[c] + int(pos) - 1
                expected = wt_seq[int(pos) - 1]
                actual = resname_at.get((c, resno))
                if actual != expected:
                    raise ValueError(
                        f"residue-identity mismatch on chain {c} position {pos} (residue {resno}): "
                        f"scores say wild-type {expected!r} but structure has {actual!r}. "
                        "The position->residue-number mapping is not valid for this structure."
                    )
                resno_value[(c, resno)] = float(value)

        # Color normalization (centered, symmetric) over all selected positions.
        all_vals = np.concatenate([s.to_numpy() for s in means.values()])
        center = 0.0 if self.ener_type == "ddG" else float(np.nanmean(all_vals))
        if vmin is None or vmax is None:
            spread = float(np.nanmax(np.abs(all_vals - center))) or 1.0
            vmin = center - spread if vmin is None else vmin
            vmax = center + spread if vmax is None else vmax

        # Sample the matplotlib colormap into an NGL colorScale array. We color via the structure's
        # B-factor + NGL's built-in "bfactor" scheme rather than a ColormakerRegistry selection
        # scheme, because a freshly-built widget can render before a just-registered custom scheme
        # syncs to the frontend (NGL then falls back to default per-chain colors).
        if cmap is None:
            # Blue -> light gray -> red, matching the heatmaps' diverging scale.
            cmap = mcolors.LinearSegmentedColormap.from_list("ddg_bgr", ["blue", "#e6e6e6", "red"])
        elif isinstance(cmap, str):
            import matplotlib

            cmap = matplotlib.colormaps[cmap]
        color_scale = [mcolors.to_hex(cmap(i / 31)) for i in range(32)]

        # Write each scored residue's value into the B-factor of all its atoms (others -> 0, but
        # they're drawn by the gray base representation so their value is irrelevant).
        b_factor = np.array(
            [resno_value.get((str(ch), int(r)), 0.0) for ch, r in zip(arr.chain_id, arr.res_id)],
            dtype=float,
        )
        arr.set_annotation("b_factor", b_factor)
        buf = io.StringIO()
        out_pdb = pdb_io.PDBFile()
        out_pdb.set_structure(arr)
        out_pdb.write(buf)

        scored_sele = " or ".join(f"{resno}:{c}" for (c, resno) in resno_value)

        view = nv.NGLWidget()
        # default_representation=False suppresses NGL's auto per-chain cartoon; otherwise it loads
        # asynchronously and clear_representations() can race it, leaving a chain-colored structure
        # rendered underneath ours.
        view.add_component(
            nv.TextStructure(buf.getvalue()), ext="pdb", default_representation=False
        )
        # Gray base for un-scored residues; colored overlay for scored residues only.
        base_sele = f"not ({scored_sele})" if scored_sele else "*"
        view.add_representation(representation, selection=base_sele, color=missing_color)
        if scored_sele:
            view.add_representation(
                representation,
                selection=scored_sele,
                color_scheme="bfactor",
                color_scale=color_scale,
                color_domain=[vmin, vmax],
            )
        view.center()
        return view


def run_energy_prediction(
    pdb_path: str | Path,
    scan_chains: list[str] | None = None,
    ddG: bool = True,
    mean_norm: bool = False,
    max_tokens: int = 20000,
    filter: bool = False,
    skip_gaps: bool = False,
    noise: float = 0.0,
    binding_partitions: list[list[str]] | None = None,
    binding_energy_cutoff: float = 8,
    mutant_csv: str | Path | None = None,
    mutant_fasta: str | Path | None = None,
    model_weights_folder: str = "vanilla_model_weights",
    model_name: str = "pottsmpnn_msa_20",
    dev: str = "cuda",
    *,
    pottsmpnn_repo: str | Path = config.POTTSMPNN_REPO,
    out_dir: str | Path | None = None,
    extra_config: Optional[dict] = None,
) -> EnergyPredictionResult:
    """Predict PottsMPNN mutation energies for a PDB and return a plottable result.

    By default (``mutant_csv`` and ``mutant_fasta`` both None) this runs a full deep
    mutational scan: every position x every non-wild-type amino acid in every chain
    (or just ``scan_chains`` if given). Provide ``mutant_csv`` or ``mutant_fasta`` to
    score a specific set of mutants instead.

    Parameters
    ----------
    pdb_path:
        Input structure. The file must be named ``<name>.pdb``; ``<name>`` becomes the
        identifier used throughout.
    scan_chains:
        Restrict the deep mutational scan to these chains (default: scan all chains). This
        only changes *which* chains are mutated and emitted — the structure is always
        encoded as the full complex, so predictions for the scanned chains are identical
        whether or not the others are scanned. Use it to avoid scoring chains you don't care
        about on large complexes.
    ddG:
        If True (default), predict relative ddG (wild type centered at 0); if False,
        predict absolute energy.
    mean_norm, max_tokens, filter, skip_gaps, noise:
        Passed through to PottsMPNN inference (see PottsMPNN docs).
    binding_partitions:
        Predict the effect of mutations on *binding* instead of stability. A list of chain
        partitions, e.g. ``[["A"], ["B"]]``, that splits the complex into its unbound groups;
        the predicted value becomes binding ddG = E(complex) - sum(E(separated partitions)).
        The partitions must cover every chain in the structure exactly once (every chain goes
        in exactly one group). Requires ``ddG=True``.
    binding_energy_cutoff:
        Angstrom cutoff used for binding-energy residues.
    mutant_csv, mutant_fasta:
        Optional explicit mutant inputs (mutually exclusive). When provided, the DMS scan
        is skipped and only these mutants are scored.
    model_weights_folder, model_name:
        Checkpoint location, relative to the PottsMPNN repo. Default
        ``vanilla_model_weights/pottsmpnn_msa_20.pt`` (vocab auto-set to 22).
    dev:
        Torch device ("cuda" or "cpu").
    pottsmpnn_repo:
        PottsMPNN repo/submodule root. Defaults to ``config.POTTSMPNN_REPO``.
    out_dir:
        If provided, PottsMPNN's output files (scores CSV, stats CSV, heatmap PNG) are
        written here and kept. If None, a temporary directory is used and discarded; the
        returned result still holds all data in memory.
    extra_config:
        Optional dict merged into the generated OmegaConf config (last), to override any
        key not exposed as an argument.

    Returns
    -------
    EnergyPredictionResult
    """
    repo_root = Path(pottsmpnn_repo).expanduser().resolve()
    if not repo_root.exists():
        raise FileNotFoundError(
            f"PottsMPNN repo not found at {repo_root}. Set POTTSMPNN_REPO or update executables.yaml."
        )

    pdb_path = Path(pdb_path).expanduser().resolve()
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    if mutant_csv is not None and mutant_fasta is not None:
        raise ValueError("Provide at most one of mutant_csv / mutant_fasta, not both.")

    pdb_name = pdb_path.stem

    # Resolve and validate the checkpoint.
    weights_dir = (repo_root / model_weights_folder).resolve()
    if not weights_dir.exists():
        raise FileNotFoundError(f"Model weights folder not found: {weights_dir}")
    checkpoint_path = (weights_dir / f"{model_name}.pt").resolve()
    if not checkpoint_path.exists():
        available = sorted(p.name for p in weights_dir.glob("*.pt"))
        hint = f"Available checkpoints: {', '.join(available)}" if available else "No .pt files found."
        raise FileNotFoundError(f"Checkpoint not found: {checkpoint_path}\n{hint}")

    # Parse the structure's chains up front (cheap) so we can validate scan_chains /
    # binding_partitions before the expensive run, and reuse the result below instead of
    # re-parsing after.
    with _temporary_sys_path(repo_root):
        from potts_mpnn_utils import parse_PDB_seq_only

        wt_info = parse_PDB_seq_only(str(pdb_path), skip_gaps=skip_gaps)
    chain_order = list(wt_info["chain_order"])
    wildtype_sequence = ":".join(wt_info[f"seq_chain_{c}"] for c in chain_order)

    # scan_chains -> PottsMPNN's exclude_chains (its complement); the DMS loop skips excluded
    # chains (run_utils.py:755), so unscanned chains are never scored or emitted.
    if scan_chains is not None:
        if not scan_chains:
            raise ValueError("scan_chains is empty: nothing would be scanned.")
        unknown = [c for c in scan_chains if c not in chain_order]
        if unknown:
            raise ValueError(f"scan_chains {unknown} not in structure chains {chain_order}.")
        exclude_chains = [c for c in chain_order if c not in scan_chains]
    else:
        exclude_chains = []

    # binding_partitions -> {pdb_name: partitions}. Binding ddG = E(complex) - sum(E(partitions)),
    # so the partitions must cover every chain exactly once (validated here, before the run, rather
    # than via PottsMPNN's deep assert at run_utils.py:777).
    if binding_partitions is not None:
        if not ddG:
            raise ValueError(
                "binding_partitions requires ddG=True (binding ddG needs the wild-type reference)."
            )
        flat = [c for part in binding_partitions for c in part]
        if sorted(flat) != sorted(chain_order):
            raise ValueError(
                "binding_partitions must cover every chain in the structure exactly once. "
                f"Structure chains: {chain_order}; partitions cover: {flat}."
            )
        binding_energy_json = {pdb_name: [list(p) for p in binding_partitions]}
    else:
        binding_energy_json = None

    input_params = {
        "pdb_path": str(pdb_path),
        "scan_chains": scan_chains,
        "ddG": ddG,
        "mean_norm": mean_norm,
        "max_tokens": max_tokens,
        "filter": filter,
        "skip_gaps": skip_gaps,
        "noise": noise,
        "binding_partitions": binding_partitions,
        "binding_energy_cutoff": binding_energy_cutoff,
        "mutant_csv": str(mutant_csv) if mutant_csv is not None else None,
        "mutant_fasta": str(mutant_fasta) if mutant_fasta is not None else None,
        "model_weights_folder": model_weights_folder,
        "model_name": model_name,
        "dev": dev,
    }

    tmp_dir = tempfile.TemporaryDirectory()
    tmp_path = Path(tmp_dir.name)

    input_list_path = tmp_path / "input_list.txt"
    input_list_path.write_text(f"{pdb_name}\n", encoding="utf-8")

    out_base = Path(out_dir).expanduser().resolve() if out_dir is not None else tmp_path
    out_base.mkdir(parents=True, exist_ok=True)

    default_model = {
        "check_path": str(checkpoint_path),
        "hidden_dim": 128,
        "edge_features": 128,
        "potts_dim": 400,
        "num_layers": 3,
        "num_edges": 48,
        "vocab": 22,  # overwritten by energy_prediction from the checkpoint name ("msa" -> 22)
    }
    default_inference = {
        "ddG": ddG,
        "mean_norm": mean_norm,
        "max_tokens": max_tokens,
        "filter": filter,
        "binding_energy_json": binding_energy_json,
        "binding_energy_cutoff": binding_energy_cutoff,
        "skip_gaps": skip_gaps,
        "noise": noise,
        "chain_dict": None,
        # chain_ranges drives the auto-plot's zoom; we keep it off here and expose zooming
        # via EnergyPredictionResult.plot_heatmap instead.
        "chain_ranges": None,
        # Empty list is a no-op in process_data's DMS loop (chain in [] -> False).
        "exclude_chains": exclude_chains or [],
    }

    cfg = OmegaConf.create(
        {
            "dev": dev,
            "out_dir": str(out_base),
            "out_name": pdb_name,
            "input_list": str(input_list_path),
            "input_dir": str(pdb_path.parent),
            "mutant_fasta": str(Path(mutant_fasta).expanduser().resolve()) if mutant_fasta is not None else None,
            "mutant_csv": str(Path(mutant_csv).expanduser().resolve()) if mutant_csv is not None else None,
            "model": default_model,
            "inference": default_inference,
        }
    )
    if extra_config:
        cfg = OmegaConf.merge(cfg, OmegaConf.create(extra_config))

    cfg_path = tmp_path / "energy_prediction_autogen.yaml"
    OmegaConf.save(cfg, cfg_path)

    energy_args = SimpleNamespace(config=str(cfg_path))

    with _temporary_sys_path(repo_root), _temporary_cwd(repo_root):
        # Imported here (not at module top) so PottsMPNN's flat bare imports resolve
        # against the repo root that _temporary_sys_path just put on sys.path. The repo's
        # module is the bare top-level ``energy_prediction``; this wrapper module is the
        # dotted ``snekwrap.wrappers.pottsmpnn.energy_prediction``, so they do not collide.
        from energy_prediction import energy_prediction as _run_energy_prediction

        _run_energy_prediction(energy_args)

    scores_path = out_base / f"{pdb_name}_scores.csv"
    if not scores_path.exists():
        raise FileNotFoundError(
            f"Expected scores file not produced by energy_prediction: {scores_path}"
        )
    scores = pd.read_csv(scores_path)
    missing = [c for c in _SCORE_COLUMNS if c not in scores.columns]
    if missing:
        raise ValueError(
            f"Scores file {scores_path} is missing expected columns {missing}; got {list(scores.columns)}"
        )

    stats_path = out_base / f"{pdb_name}_stats.csv"
    stats = pd.read_csv(stats_path) if stats_path.exists() else None

    return EnergyPredictionResult(
        pdb_name=pdb_name,
        scores=scores,
        chain_order=chain_order,
        wildtype_sequence=wildtype_sequence,
        ener_type="ddG" if ddG else "dG",
        repo_root=repo_root,
        input_params=input_params,
        stats=stats,
    )
