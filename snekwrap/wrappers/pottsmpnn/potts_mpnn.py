"""
Thin functional wrapper around PottsMPNN so it can be called from Python
(similar to ``run_protein_mpnn``).

The function expects a PottsMPNN OmegaConf YAML config (same format used by
``sample_seqs.py``) and returns parsed sequences and metrics as Python objects
instead of writing a CLI.
"""
# %%
from __future__ import annotations

import contextlib
import csv
import json
import os
import sys
import tempfile
import shutil
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, Iterable, Optional, Tuple

from omegaconf import OmegaConf

import snekwrap.config as config

# %%
# NOTE: PottsMPNN's science modules (sample_seqs, run_utils, potts_mpnn_utils,
# etab_utils) are NOT vendored here. They live in the PottsMPNN repo/submodule
# (config.POTTSMPNN_REPO) and use flat, bare imports, so `sample_seqs` is imported
# lazily inside run_potts_mpnn() while that repo root is on sys.path.


@dataclass
class PottsMPNNSample:
    """Container for a PottsMPNN sequence sample and metrics."""

    sequence: str  # full joined sequence (':' between chains)
    chain_sequences: Dict[str, str]  # per-chain sequences
    sample_number: Optional[int] = None
    source_file: Optional[Path] = None
    optimization_mode: Optional[str] = None
    optimized_sequence: Optional[str] = None
    optimized_chain_sequences: Optional[Dict[str, str]] = None
    optimized_source_file: Optional[Path] = None
    # PottsMPNN whole-sequence Potts energy of `sequence` (arbitrary units), only set when
    # run_potts_mpnn(..., compute_energies=True). See _compute_potts_energies for caveats.
    energy: Optional[float] = None
    # Energy of `optimized_sequence` (a.u.); set only if optimization ran and compute_energies=True.
    optimized_energy: Optional[float] = None


@contextlib.contextmanager
def _temporary_sys_path(path: Path) -> Iterable[None]:
    """Temporarily prepend a path to sys.path while the context is active."""

    path_str = str(path)
    sys.path.insert(0, path_str)
    try:
        yield
    finally:
        try:
            sys.path.remove(path_str)
        except ValueError:
            pass


@contextlib.contextmanager
def _temporary_cwd(path: Path) -> Iterable[None]:
    """Temporarily switch the working directory."""

    original = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(original)


def _resolve_repo_path(raw_path: str | Path, repo_root: Path) -> Path:
    """Resolve a path, treating relative inputs as rooted at the repo."""

    resolved = Path(raw_path)
    if not resolved.is_absolute():
        resolved = (repo_root / resolved).resolve()
    return resolved


def _read_fasta(path: Path) -> Dict[str, str]:
    sequences: Dict[str, str] = {}
    with path.open() as handle:
        lines = [line.strip() for line in handle if line.strip()]
    for i in range(0, len(lines), 2):
        header = lines[i].lstrip(">")
        if i + 1 >= len(lines):
            break
        sequences[header] = lines[i + 1]
    return sequences


def _read_metrics(path: Path) -> Dict[str, Dict[str, float]]:
    metrics: Dict[str, Dict[str, float]] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            name = row.get("pdb") or row.get("name") or ""
            if not name:
                continue
            metrics[name] = {}
            for key in ("seq_loss", "nsr", "potts_loss"):
                value = row.get(key)
                if value is None or value == "":
                    metrics[name][key] = None
                else:
                    try:
                        metrics[name][key] = float(value)
                    except ValueError:
                        metrics[name][key] = None
    return metrics


def _parse_sample_number(name: str) -> Optional[int]:
    suffix = name.split("_")[-1]
    return int(suffix) if suffix.isdigit() else None


def _split_chain_sequences(seq: str, chain_order: list[str]) -> Dict[str, str]:
    parts = seq.split(":")
    if len(parts) == len(chain_order):
        return dict(zip(chain_order, parts))

    # Fallback: still return what we have with best-effort keys
    chain_sequences: Dict[str, str] = {}
    for idx, part in enumerate(parts):
        key = chain_order[idx] if idx < len(chain_order) else f"chain{idx + 1}"
        chain_sequences[key] = part
    return chain_sequences


def _cleanup_outputs(out_dir: Path, out_name: str, optimization_mode: str | None = None) -> None:
    """Remove on-disk artifacts produced by PottsMPNN for a single run."""

    files = [
        out_dir / f"{out_name}.fasta",
        out_dir / f"{out_name}_decoding_order.json",
        out_dir / f"{out_name}_av_loss.csv",
    ]

    if optimization_mode:
        files.append(out_dir / f"{out_name}_optimized_{optimization_mode}.fasta")

    pdb_dir = out_dir / f"{out_name}_pdbs"

    for path in files:
        try:
            path.unlink(missing_ok=True)
        except OSError:
            pass

    if pdb_dir.exists():
        shutil.rmtree(pdb_dir, ignore_errors=True)


def _compute_potts_energies(
    cfg,
    repo_root: Path,
    chain_list: list[str],
    sequences: list[str],
) -> list[float]:
    """Recompute PottsMPNN whole-sequence energies for already-generated sequences.

    ``sample_seqs`` ranks its samples by the Potts energy
    ``E(seq) = sum_i h_i(s_i) + sum_{i<j} J_ij(s_i, s_j)`` but never writes it out. This
    reloads the model, encodes the backbone once into a single energy table (etab), and
    scores each sequence against it via ``calc_eners`` -- the same computation
    ``sample_seqs`` performs internally.

    The returned energies are in arbitrary units. They are comparable *across these
    sequences* because they all share one backbone/etab (lower = the model prefers that
    sequence on this structure), but they are NOT comparable across different structures.

    ``sequences`` are the ``':'``-joined per-chain sequences with chains in ``chain_list``
    order; the ``':'`` separators are stripped before scoring so positions align 1:1 with
    the etab.
    """
    with _temporary_sys_path(repo_root), _temporary_cwd(repo_root):
        # Bare imports resolve against the repo root just added to sys.path.
        import torch
        from run_utils import get_etab
        from potts_mpnn_utils import PottsMPNN, parse_PDB
        import etab_utils

        # sample_seqs derives the vocab size from the checkpoint name (msa models use 22).
        # Replicate it so the model and featurization match the run that produced these
        # sequences; the wrapper's in-memory cfg.model.vocab is not authoritative.
        vocab = 22 if "msa" in cfg.model.check_path else 21

        checkpoint = torch.load(cfg.model.check_path, map_location="cpu", weights_only=False)
        model = PottsMPNN(
            ca_only=False,
            num_letters=vocab,
            vocab=vocab,
            node_features=cfg.model.hidden_dim,
            edge_features=cfg.model.hidden_dim,
            hidden_dim=cfg.model.hidden_dim,
            potts_dim=cfg.model.potts_dim,
            num_encoder_layers=cfg.model.num_layers,
            num_decoder_layers=cfg.model.num_layers,
            k_neighbors=cfg.model.num_edges,
            augment_eps=cfg.inference.noise,
        )
        model.load_state_dict(checkpoint["model_state_dict"], strict=False)
        model.eval()
        model = model.to(cfg.dev)
        for param in model.parameters():
            param.requires_grad = False

        # get_etab reads cfg.model.vocab during featurization; keep it consistent.
        cfg.model.vocab = vocab

        input_pdb = os.path.join(cfg.input_dir, cfg.out_name + ".pdb")
        pdb_data = parse_PDB(input_pdb, chain_list, skip_gaps=cfg.inference.skip_gaps)
        etab, E_idx, wt_seq = get_etab(model, pdb_data, cfg, None)

        energies: list[float] = []
        for seq in sequences:
            seq_flat = seq.replace(":", "")
            assert len(seq_flat) == len(wt_seq), (
                f"generated sequence length {len(seq_flat)} != structure length "
                f"{len(wt_seq)}; cannot align sequence to the energy table"
            )
            seq_tensor = etab_utils.seq_to_tensor(seq_flat).to(
                dtype=torch.int64, device=E_idx.device
            )
            # calc_eners expects [B, n, L]; score one sequence at a time to keep the
            # expanded etab small for large complexes.
            batch_scores, _, _ = etab_utils.calc_eners(
                etab, E_idx, seq_tensor.view(1, 1, -1), None
            )
            energies.append(batch_scores.squeeze().cpu().item())
    return energies


def run_potts_mpnn(
    pdb_path: str | Path,
    designed_chains: list[str] | None = None,
    fixed_chains: list[str] | None = None,
    fixed_position_chain: str | None = None,
    fixed_positions: list[int] | None = None,
    num_seqs: int = 8,
    sampling_temp: str = "0.1",
    optimization_temperature: str = "0.0",
    model_weights_folder: str = "soluble_model_weights",
    backbone_noise: float = 0.0,
    model_name: str = "sol_pottsmpnn_msa_20",
    omit_AAs: list[str] | None = None,
    dev: str = "cuda",
    compute_energies: bool = False,
    *,
    pottsmpnn_repo: str | Path = config.POTTSMPNN_REPO,
    out_dir: str | Path | None = None,
    extra_config: Optional[dict] = None,
) -> Tuple[list[PottsMPNNSample], list[PottsMPNNSample], Dict]:
    """
    Execute PottsMPNN sequence generation with a call signature mirroring
    ``run_protein_mpnn`` (single PDB path, chain selections, sampling params).

    Returns a tuple of (samples, optimized_samples, input_params). The first element always
    contains the raw sampled sequences. If optimization_mode is set, the second
    element contains optimized sequences (also attached to the matching sample
    objects when names align). The third element is a dictionary of all input parameters.

    Parameters
    ----------
    pdb_path: PDB file to design.
    designed_chains: Chains to design (defaults to ["A"]).
    fixed_chains: Chains to keep fixed (defaults to []).
    fixed_position_chain: Chain on which specific positions are fixed.
    fixed_positions: 1-based residue indices to fix within fixed_position_chain.
    num_seqs: Number of sequences to sample.
    sampling_temp: Sampling temperature string (space-separated allowed).
    model_weights_folder: Folder under the PottsMPNN repo containing weights.
    backbone_noise: Backbone noise (maps to PottsMPNN "noise").
    model_name: Checkpoint basename (without extension).
    compute_energies: If True, recompute each sampled (and optimized) sequence's PottsMPNN
        whole-sequence energy and attach it to the sample's ``energy`` (and
        ``optimized_energy``) field. Off by default because it reloads the model and runs an
        extra encoder pass. Energies are arbitrary units, comparable only across sequences
        on this same backbone (see ``_compute_potts_energies``).
    pottsmpnn_repo: Location of the PottsMPNN repository.
    out_dir: Optional override for output directory (default: outputs/pottsmpnn).
    extra_config: Optional dict merged into the generated config.
    """

    # Capture input parameters
    input_params = {
        "pdb_path": str(pdb_path),
        "designed_chains": designed_chains or ["A"],
        "fixed_chains": fixed_chains or [],
        "fixed_position_chain": fixed_position_chain,
        "fixed_positions": fixed_positions,
        "num_seqs": num_seqs,
        "sampling_temp": sampling_temp,
        "optimization_temperature": optimization_temperature,
        "model_weights_folder": model_weights_folder,
        "backbone_noise": backbone_noise,
        "model_name": model_name,
        "omit_AAs": omit_AAs or [],
        "dev": dev,
        "compute_energies": compute_energies,
        "pottsmpnn_repo": str(pottsmpnn_repo),
        "out_dir": str(out_dir or "outputs/pottsmpnn"),
        "extra_config": extra_config,
    }

    repo_root = Path(pottsmpnn_repo).expanduser().resolve()
    if not repo_root.exists():
        raise FileNotFoundError(
            f"PottsMPNN repo not found at {repo_root}. Set POTTSMPNN_REPO or update executables.yaml."
        )

    pdb_path = Path(pdb_path).expanduser().resolve()
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    designed_chains = designed_chains or ["A"]
    fixed_chains = fixed_chains or []
    chain_list = list(set(designed_chains + fixed_chains))
    chain_order = designed_chains + fixed_chains  # preserve caller-provided ordering

    # Build chain dict JSON (PottsMPNN expects a mapping of pdb -> [design, fixed])
    tmp_dir = tempfile.TemporaryDirectory()
    tmp_path = Path(tmp_dir.name)

    pdb_name = pdb_path.stem
    chain_dict = {pdb_name: [designed_chains, fixed_chains]}
    chain_dict_json = tmp_path / "chain_dict.json"
    # PottsMPNN expects a JSON file; write explicitly to avoid YAML formatting
    chain_dict_json.write_text(json.dumps(chain_dict), encoding="utf-8")

    # Fixed positions dict if provided
    fixed_positions_json = ""
    if fixed_positions:
        if fixed_position_chain is None:
            raise ValueError("fixed_position_chain must be provided when fixed_positions are set")
        fixed_positions_dict = {pdb_name: {c: [] for c in chain_list}}
        fixed_positions_dict[pdb_name][fixed_position_chain] = [int(p) for p in fixed_positions]
        fixed_positions_json = tmp_path / "fixed_positions.json"
        # process_configs reads JSON lines; emit compact JSON on a single line
        fixed_positions_json.write_text(json.dumps(fixed_positions_dict), encoding="utf-8")

    # Input list file (one PDB entry)
    input_list_path = tmp_path / "input_list.txt"
    input_dir = pdb_path.parent
    input_list_path.write_text(f"{pdb_name}\n", encoding="utf-8")

    # Model/inference defaults
    weights_dir = (repo_root / model_weights_folder).resolve()
    if not weights_dir.exists():
        raise FileNotFoundError(f"Model weights folder not found: {weights_dir}")

    checkpoint_path = (weights_dir / f"{model_name}.pt").resolve()
    if not checkpoint_path.exists():
        available = sorted(p.name for p in weights_dir.glob("*.pt"))
        hint = f"Available checkpoints: {', '.join(available)}" if available else "No .pt files found in weights folder."
        raise FileNotFoundError(f"Checkpoint not found: {checkpoint_path}\n{hint}")

    if omit_AAs is None:
        omit_AAs = []
    default_model = {
        "check_path": str(checkpoint_path),
        "hidden_dim": 128,
        "edge_features": 128,
        "potts_dim": 400,
        "num_layers": 3,
        "num_edges": 48,
        "vocab": 21,
    }
    default_inference = {
        "num_samples": num_seqs,
        "temperature": float(str(sampling_temp).split()[0]),
        # Match sampling temp unless caller overrides; required by sample_seqs
        "optimization_temperature": float(str(optimization_temperature).split()[0]),
        "noise": backbone_noise,
        "skip_gaps": False,
        "fix_decoding_order": True,
        "decoding_order_offset": 0,
        "optimization_mode": "potts",
        "binding_energy_optimization": "none",
        "binding_energy_json": None,
        "binding_energy_cutoff": 8,
        "optimize_pdb": False,
        "optimize_fasta": "",
        # PDBs written by sample_seqs are not surfaced by this wrapper (they're
        # cleaned up unused), so skip writing them.
        "write_pdb": False,
        "fixed_positions_json": str(fixed_positions_json),
        "pssm_json": "",
        "omit_AA_json": "",
        "bias_AA_json": "",
        "tied_positions_json": "",
        "tied_epistasis": True,
        "bias_by_res_json": "",
        "omit_AAs": omit_AAs,
        "pssm_threshold": 0.0,
        "pssm_multi": 0.0,
        "pssm_log_odds_flag": False,
        "pssm_bias_flag": False,
    }

    cfg = OmegaConf.create(
        {
            "dev": dev,
            "out_dir": str(out_dir or "outputs/pottsmpnn"),
            "out_name": pdb_name,
            "input_list": str(input_list_path),
            "input_dir": str(input_dir),
            "chain_dict_json": str(chain_dict_json),
            "model": default_model,
            "inference": default_inference,
        }
    )

    if extra_config:
        cfg = OmegaConf.merge(cfg, OmegaConf.create(extra_config))

    cfg.model.check_path = str(_resolve_repo_path(cfg.model.check_path, repo_root))

    final_checkpoint = Path(cfg.model.check_path)
    if not final_checkpoint.exists():
        raise FileNotFoundError(f"Resolved checkpoint not found: {final_checkpoint}")

    cfg_path = tmp_path / "pottsmpnn_autogen.yaml"
    OmegaConf.save(cfg, cfg_path)

    out_dir_path = _resolve_repo_path(cfg.out_dir, repo_root)
    out_dir_path.mkdir(parents=True, exist_ok=True)

    sample_args = SimpleNamespace(config=str(cfg_path))

    with _temporary_sys_path(repo_root), _temporary_cwd(repo_root):
        # Imported here (not at module top) so PottsMPNN's flat bare imports
        # (`from run_utils import ...`, `import etab_utils`) resolve against the
        # repo root that _temporary_sys_path just put on sys.path.
        from sample_seqs import sample_seqs

        sample_seqs(sample_args)

    samples: list[PottsMPNNSample] = []
    samples_by_name: Dict[str, PottsMPNNSample] = {}
    optimized_samples: list[PottsMPNNSample] = []

    fasta_path = out_dir_path / f"{cfg.out_name}.fasta"
    metrics_path = out_dir_path / f"{cfg.out_name}_av_loss.csv"
    metrics = _read_metrics(metrics_path) if metrics_path.exists() else {}

    if fasta_path.exists():
        for name, seq in _read_fasta(fasta_path).items():
            chain_sequences = _split_chain_sequences(seq, chain_order)

            sample = PottsMPNNSample(
                sequence=seq,
                chain_sequences=chain_sequences,
                sample_number=_parse_sample_number(name),
                source_file=fasta_path,
            )
            samples.append(sample)
            samples_by_name[name] = sample

    if getattr(cfg, "inference", None) and getattr(cfg.inference, "optimization_mode", ""):
        opt_path = out_dir_path / f"{cfg.out_name}_optimized_{cfg.inference.optimization_mode}.fasta"
        if opt_path.exists():
            for name, seq in _read_fasta(opt_path).items():
                chain_sequences = _split_chain_sequences(seq, chain_order)

                optimization_mode = getattr(cfg.inference, "optimization_mode", None) or None
                optimized_sample = PottsMPNNSample(
                    sequence=seq,
                    chain_sequences=chain_sequences,
                    sample_number=_parse_sample_number(name),
                    source_file=opt_path,
                    optimization_mode=optimization_mode,
                    optimized_sequence=seq,
                    optimized_chain_sequences=chain_sequences,
                    optimized_source_file=opt_path,
                )

                optimized_samples.append(optimized_sample)

                if name in samples_by_name:
                    # Also attach to the originating sample for convenience
                    existing = samples_by_name[name]
                    existing.optimization_mode = optimization_mode
                    existing.optimized_sequence = seq
                    existing.optimized_chain_sequences = chain_sequences
                    existing.optimized_source_file = opt_path

    if compute_energies:
        # Score every sequence against a single freshly-built energy table. Gather all the
        # (object, attribute, sequence) targets first so the model is loaded / encoded once.
        energy_targets: list[tuple[PottsMPNNSample, str, str]] = []
        for sample in samples:
            energy_targets.append((sample, "energy", sample.sequence))
            if sample.optimized_sequence is not None:
                energy_targets.append((sample, "optimized_energy", sample.optimized_sequence))
        for optimized_sample in optimized_samples:
            energy_targets.append((optimized_sample, "energy", optimized_sample.sequence))

        if energy_targets:
            energies = _compute_potts_energies(
                cfg, repo_root, chain_order, [seq for _, _, seq in energy_targets]
            )
            for (obj, attr, _), energy in zip(energy_targets, energies, strict=True):
                setattr(obj, attr, energy)

    _cleanup_outputs(out_dir_path, cfg.out_name, getattr(cfg.inference, "optimization_mode", None))
    tmp_dir.cleanup()
    return samples, optimized_samples, input_params
