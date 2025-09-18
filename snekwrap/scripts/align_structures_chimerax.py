#!/Applications/ChimeraX-1.9.app/Contents/bin/ChimeraX --nogui


# %%
import os
from pathlib import Path
# %%

# PDB_DIRECTORY = Path(__file__).resolve().parent

from chimerax.core.commands import run
from chimerax.atomic import all_atomic_structures
import argparse
PDB_DIRECTORY = Path.cwd()


def generate_sessions(pdb_files, filename):
    # Open the reference PDB file (assumed to be the first one)
    run(session, f"open '{pdb_files[0]}'")
    reference_model_id = 1  # The first model (assumed as #1)
    model_num_map = {reference_model_id: pdb_files[0]}  # Map model ID to PDB file path
    # Align each subsequent PDB file to the reference (first) model
    for i, pdb_file in enumerate(pdb_files[1:], start=2):
        run(session, f"open '{pdb_file}'")
        # Align the newly opened model (model id will be `i`) to the reference model (model id #1)
        run(session, f"mm #{i} to #{reference_model_id}")
        model_num_map[i] = pdb_file  # Map model ID to PDB file path
    run(session, "color bychain cartoons")
    # run(session, "rainbow ~/A residues")
    run(session, "transparency 50 cartoons")
    run(session, "transparency 50 atoms")
    run(session, "color /A gray")
    # run(session, "color bfactor /B range 65,100 palette blue:red; key blue:65 red:100 size 0.3,0.05 pos 0.4,0.1")
    # zoom out a bit
    run(session, "view")
    run(session, "set bgColor white")
    run(session, "graphics silhouettes true")
    run(session, f"save '{filename}'")
    run(session, f"close all")

def main(pdb_directory):
    pdb_files = list(pdb_directory.glob('*.pdb'))
    filename = pdb_directory.resolve().name + ".cxs"
    generate_sessions(pdb_files, filename)
    run(session, "exit")

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="Align structures in ChimeraX. saves as .cxs file", formatter_class=argparse.RawDescriptionHelpFormatter)
#     parser.add_argument(
#         "--pdb_directory",
#         type=str,
#         help="Directory containing PDB files. This is also where the .cxs file will be saved. Default is current directory",
#         metavar="<file>",
#         default=PDB_DIRECTORY,
#     )
#     args = parser.parse_args()
#     main(args.pdb_directory)

main(PDB_DIRECTORY)