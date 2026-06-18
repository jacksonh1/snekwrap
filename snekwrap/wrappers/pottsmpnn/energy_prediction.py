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
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List, Optional

import pandas as pd
from omegaconf import OmegaConf

import snekwrap.config as config

# Reuse the sys.path / cwd helpers from the sequence-design wrapper rather than
# duplicating them. Importing potts_mpnn does not import this module, so there is no
# circular import.
from snekwrap.wrappers.pottsmpnn.potts_mpnn import _temporary_sys_path, _temporary_cwd


_SCORE_COLUMNS = ["pdb", "mutant", "wildtype", "ddG_pred", "ddG_expt"]


@dataclass
class EnergyPredictionResult:
    """PottsMPNN energy-prediction results for a single PDB.

    The mutation scores live in ``scores`` (one row per scored mutant). Use
    :meth:`plot_heatmap` to render the per-position x per-amino-acid heatmap.
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
            subset of a large complex is scanned (e.g. ``exclude_chains``), to avoid hundreds
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


def run_energy_prediction(
    pdb_path: str | Path,
    exclude_chains: list[str] | None = None,
    ddG: bool = True,
    mean_norm: bool = False,
    max_tokens: int = 20000,
    filter: bool = False,
    skip_gaps: bool = False,
    noise: float = 0.0,
    binding_energy_json: dict | str | None = None,
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
    (minus ``exclude_chains``). Provide ``mutant_csv`` or ``mutant_fasta`` to score a
    specific set of mutants instead.

    Parameters
    ----------
    pdb_path:
        Input structure. The file must be named ``<name>.pdb``; ``<name>`` becomes the
        identifier used throughout.
    exclude_chains:
        Chains to skip during the deep mutational scan.
    ddG:
        If True (default), predict relative ddG (wild type centered at 0); if False,
        predict absolute energy.
    mean_norm, max_tokens, filter, skip_gaps, noise:
        Passed through to PottsMPNN inference (see PottsMPNN docs).
    binding_energy_json:
        Optional chain-partition spec (dict or path to JSON) to predict the effect of
        mutations on binding energy instead of stability.
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

    input_params = {
        "pdb_path": str(pdb_path),
        "exclude_chains": exclude_chains or [],
        "ddG": ddG,
        "mean_norm": mean_norm,
        "max_tokens": max_tokens,
        "filter": filter,
        "skip_gaps": skip_gaps,
        "noise": noise,
        "binding_energy_json": binding_energy_json,
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
        from potts_mpnn_utils import parse_PDB_seq_only

        _run_energy_prediction(energy_args)

        wt_info = parse_PDB_seq_only(str(pdb_path), skip_gaps=skip_gaps)

    chain_order = list(wt_info["chain_order"])
    wildtype_sequence = ":".join(wt_info[f"seq_chain_{c}"] for c in chain_order)

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
