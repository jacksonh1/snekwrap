#!/Applications/ChimeraX-1.9.app/Contents/bin/ChimeraX --nogui

import os
from chimerax.core.commands import run
from chimerax.atomic import all_atomic_structures
from pathlib import Path
PDB_DIRECTORY = Path.cwd()



output_dir = PDB_DIRECTORY / 'aligned_pdbs'
output_dir.mkdir(exist_ok=True, parents=True)
pdb_files = list(PDB_DIRECTORY.glob('*.pdb'))
filename = PDB_DIRECTORY.resolve().name + ".cxs"

run(session, f"open '{pdb_files[0]}'")
reference_model_id = 1  # The first model (assumed as #1)
model_num_map = {reference_model_id: pdb_files[0]}  # Map model ID to PDB file path

# Align each subsequent PDB file to the reference (first) model
for i, pdb_file in enumerate(pdb_files[1:], start=2):
    run(session, f"open '{pdb_file}'")
    # Align the newly opened model (model id will be `i`) to the reference model (model id #1)
    run(session, f"mm #{i} to #{reference_model_id}")
    model_num_map[i] = pdb_file  # Map model ID to PDB file path

# Save aligned models
for model in all_atomic_structures(session):
    model_id = model.id[0]
    print(f"Saving model with ID: {model_id}")
    run(session, f"save '{str(output_dir)}/{model_num_map[model_id].stem}-aligned.pdb' format pdb models #{model_id}")

run(session, f"save '{filename}'")
run(session, "exit")

