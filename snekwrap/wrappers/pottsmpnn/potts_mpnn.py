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
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, Iterable, Optional

from omegaconf import OmegaConf

import snekwrap.config as config

# %%

from snekwrap.wrappers.pottsmpnn.sample_seqs import sample_seqs


@dataclass
class PottsMPNNSample:
    """Container for a PottsMPNN sequence sample and metrics."""

    sequence: str  # full joined sequence (':' between chains)
    chain_sequences: Dict[str, str]  # per-chain sequences
    sample_number: Optional[int] = None
    source_file: Optional[Path] = None


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


def run_potts_mpnn(
    pdb_path: str | Path,
    designed_chains: list[str] | None = None,
    fixed_chains: list[str] | None = None,
    fixed_position_chain: str | None = None,
    fixed_positions: list[int] | None = None,
    num_seqs: int = 8,
    sampling_temp: str = "0.1",
    model_weights_folder: str = "vanilla_model_weights",
    backbone_noise: float = 0.0,
    model_name: str = "pottsmpnn_20",
    *,
    pottsmpnn_repo: str | Path = config.POTTSMPNN_REPO,
    out_dir: str | Path | None = None,
    extra_config: Optional[dict] = None,
) -> list[PottsMPNNSample]:
    """
    Execute PottsMPNN sequence generation with a call signature mirroring
    ``run_protein_mpnn`` (single PDB path, chain selections, sampling params).

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
    pottsmpnn_repo: Location of the PottsMPNN repository.
    out_dir: Optional override for output directory (default: outputs/pottsmpnn).
    extra_config: Optional dict merged into the generated config.
    """

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
        "write_pdb": True,
        "fixed_positions_json": str(fixed_positions_json),
        "pssm_json": "",
        "omit_AA_json": "",
        "bias_AA_json": "",
        "tied_positions_json": "",
        "bias_by_res_json": "",
        "omit_AAs": [],
        "pssm_threshold": 0.0,
        "pssm_multi": 0.0,
        "pssm_log_odds_flag": False,
        "pssm_bias_flag": False,
    }

    cfg = OmegaConf.create(
        {
            "dev": "cuda",
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

        sample_seqs(sample_args)

    results: list[PottsMPNNSample] = []

    fasta_path = out_dir_path / f"{cfg.out_name}.fasta"
    metrics_path = out_dir_path / f"{cfg.out_name}_av_loss.csv"
    metrics = _read_metrics(metrics_path) if metrics_path.exists() else {}

    if fasta_path.exists():
        for name, seq in _read_fasta(fasta_path).items():
            chain_sequences = _split_chain_sequences(seq, chain_order)

            results.append(
                PottsMPNNSample(
                    sequence=seq,
                    chain_sequences=chain_sequences,
                    sample_number=_parse_sample_number(name),
                    source_file=fasta_path,
                )
            )

    if getattr(cfg, "inference", None) and getattr(cfg.inference, "optimization_mode", ""):
        opt_path = out_dir_path / f"{cfg.out_name}_optimized_{cfg.inference.optimization_mode}.fasta"
        if opt_path.exists():
            for name, seq in _read_fasta(opt_path).items():
                chain_sequences = _split_chain_sequences(seq, chain_order)

                results.append(
                    PottsMPNNSample(
                        sequence=seq,
                        chain_sequences=chain_sequences,
                        sample_number=_parse_sample_number(name),
                        source_file=opt_path,
                    )
                )

    tmp_dir.cleanup()
    return results
