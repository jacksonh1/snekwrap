# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: snekwrap
#     language: python
#     name: python3
# ---

# %% [markdown]
# # PottsMPNN sequence generation
#
# `run_potts_mpnn` runs PottsMPNN's sampler over a backbone and returns designed sequences as
# Python objects (no CLI, no files to manage). It's the design counterpart to
# `run_energy_prediction`: instead of scoring mutations on a fixed sequence, it *generates* new
# sequences for the structure.
#
# Each call returns a 3-tuple `(samples, optimized_samples, input_params)`. Every design comes in
# **two flavors**:
#
# - `sample.sequence` — the raw **autoregressive** sample (ProteinMPNN-style: residues decoded one
#   at a time, each conditioned on the ones already placed).
# - `sample.optimized_sequence` — that sample after **Potts optimization**: coordinate descent on
#   the explicit energy table `E(seq) = Σ h_i + Σ J_ij`, which revisits positions to resolve
#   conflicts a single autoregressive pass can't (optimization is on by default).
#
# > PottsMPNN is bundled as the `external/PottsMPNN` submodule; no path configuration is needed.

# %%
from pathlib import Path

import snekwrap.config as config
from snekwrap.wrappers.pottsmpnn import run_potts_mpnn

# Example PDBs ship with the PottsMPNN submodule.
example_pdbs = Path(config.POTTSMPNN_REPO) / "inputs" / "example_pdbs"

# Everything here runs on CPU so the example is portable; use dev="cuda" if a GPU is available.
dev = "cpu"
print(f"running on: {dev}")

# %% [markdown]
# ## Basic generation
#
# Pass a PDB and the chain(s) to design. The default model is
# `soluble_model_weights/sol_pottsmpnn_msa_20.pt`. `3dkm` is a small single-chain protein, fast on
# CPU.

# %%
samples, optimized, params = run_potts_mpnn(
    example_pdbs / "3dkm.pdb",
    designed_chains=["A"],
    num_seqs=6,
    dev=dev,
)

print("designs generated:", len(samples))
first = samples[0]
print("sample number:", first.sample_number)
print("per-chain sequences:", {c: len(s) for c, s in first.chain_sequences.items()}, "residues")
print("sequence:", first.sequence)

# %% [markdown]
# ## Sampled vs. Potts-optimized
#
# Every sample carries both its raw autoregressive `sequence` and its `optimized_sequence`. The raw
# sample is a single left-to-right pass, so an early residue can end up in conflict with a choice
# made later; the Potts optimizer sweeps over positions and minimizes the energy table, which
# typically changes a handful of residues. The optimized set is also returned separately as
# `optimized_samples`.

# %%
raw = first.sequence
opt = first.optimized_sequence
n_changed = sum(a != b for a, b in zip(raw, opt))
print("raw sample :", raw)
print("optimized  :", opt)
print(f"positions changed by Potts optimization: {n_changed} / {len(raw.replace(':', ''))}")
print("optimized_samples returned separately:", len(optimized))

# %% [markdown]
# ## Ranking designs by predicted energy
#
# By default the wrapper doesn't score the designs (it would cost an extra encoder pass). Pass
# `compute_energies=True` to attach the PottsMPNN whole-sequence energy to each sample: `.energy`
# for the raw sample and `.optimized_energy` for its optimized version. When energies are computed,
# `samples` come back **sorted best-first** (lowest energy).
#
# These energies are in arbitrary units and are only comparable *across sequences on this same
# backbone* — lower means the model prefers that sequence on this structure. They are **not**
# comparable across different structures.

# %%
scored, _, _ = run_potts_mpnn(
    example_pdbs / "3dkm.pdb",
    designed_chains=["A"],
    num_seqs=6,
    compute_energies=True,
    dev=dev,
)

print(f"{'sample':>6} {'energy':>12} {'opt_energy':>12}")
for s in scored:
    print(f"{s.sample_number:>6} {s.energy:>12.3f} {s.optimized_energy:>12.3f}")

print("\nreturned best-first?", [round(s.energy, 2) for s in scored])
print("optimized <= raw for every sample?",
      all(s.optimized_energy <= s.energy + 1e-6 for s in scored))

# %% [markdown]
# ## Designing a subset of chains
#
# For a complex, `designed_chains` are redesigned while `fixed_chains` keep their native sequence.
# `6w25` is a peptide (chain A) bound to a receptor (chain B); here we design the peptide and hold
# the receptor fixed. Because only chain A is designed, chain B is identical across every sample.

# %%
complex_samples, _, _ = run_potts_mpnn(
    example_pdbs / "6w25.pdb",
    designed_chains=["A"],  # design the peptide
    fixed_chains=["B"],     # keep the receptor native
    num_seqs=4,
    compute_energies=True,
    dev=dev,
)

chain_a_designs = {s.chain_sequences["A"] for s in complex_samples}
chain_b_seqs = {s.chain_sequences["B"] for s in complex_samples}
print("chains per design:", list(complex_samples[0].chain_sequences))
print(f"distinct chain-A designs: {len(chain_a_designs)} / {len(complex_samples)}")
print(f"distinct chain-B sequences: {len(chain_b_seqs)} (fixed -> should be 1)")

best = complex_samples[0]  # energy-sorted, so index 0 is the best
print("best peptide design:", best.chain_sequences["A"], "| energy:", round(best.energy, 3))

# %% [markdown]
# ## Fixing specific positions
#
# To hold individual residues while designing the rest of a chain, pass `fixed_position_chain` and
# `fixed_positions` (1-based indices within that chain). Those positions keep their wild-type
# identity in every sample.

# %%
fixed_positions = [5, 10, 15]
fixed_samples, _, _ = run_potts_mpnn(
    example_pdbs / "3dkm.pdb",
    designed_chains=["A"],
    fixed_position_chain="A",
    fixed_positions=fixed_positions,
    num_seqs=4,
    dev=dev,
)

for pos in fixed_positions:
    residues = {s.chain_sequences["A"][pos - 1] for s in fixed_samples}
    print(f"position {pos}: {residues}  (fixed -> single residue across all samples)")

# %% [markdown]
# ## Controlling generation
#
# A few knobs shape the output:
#
# - **`sampling_temp`** — temperature of the autoregressive sampler. Higher values give more diverse
#   (less conservative) samples.
# - **`optimization_temperature`** — temperature of the Potts refinement step. `0` (the default) is
#   near-greedy energy minimization; higher values make the refinement stochastic.
# - **`backbone_noise`** — Gaussian noise added to the backbone coordinates in the encoder, which
#   perturbs the energy table itself (useful for sampling around an imperfect backbone).
# - **`omit_AAs`** — residue types to exclude entirely. The exclusion is enforced in *both* the
#   autoregressive sample and the Potts-optimized sequence.
# - **`num_seqs`** — how many designs to generate.
#
# Below we forbid cysteine and confirm it appears in neither the sampled nor the optimized design.

# %%
no_cys, _, _ = run_potts_mpnn(
    example_pdbs / "3dkm.pdb",
    designed_chains=["A"],
    num_seqs=4,
    omit_AAs=["C"],
    dev=dev,
)

sampled_has_cys = any("C" in s.chain_sequences["A"] for s in no_cys)
optimized_has_cys = any("C" in s.optimized_chain_sequences["A"] for s in no_cys)
print("cysteine in any sampled design?  ", sampled_has_cys)
print("cysteine in any optimized design?", optimized_has_cys)
assert not sampled_has_cys and not optimized_has_cys, "omit_AAs was not honored"
print("omit_AAs=['C'] honored in both sampled and optimized designs")
