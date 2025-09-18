import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from importlib_resources import files

ALIGN_SCRIPT_FILE = files("snekwrap.scripts").joinpath("align_structures_chimerax.py")
ALIGN_SCRIPT_FILE = Path(ALIGN_SCRIPT_FILE)
# I could copy the script into the directory and run it, or I could change
# into the directory, have the script in path, and run it.
ALIGN_EXPORT_SCRIPT_FILE = files("snekwrap.scripts").joinpath("align_and_export_pdbs_chimerax.py")
ALIGN_EXPORT_SCRIPT_FILE = Path(ALIGN_EXPORT_SCRIPT_FILE)


import snekwrap.config as config
import shutil


def align_pdbs_chimerax(
    pdb_dir: str | Path,
    chimera_executable: str = config.CHIMERAX_EXECUTABLE,
):
    """
    Changes the working directory to the pdb_dir, runs the script in there
    """
    cwd = Path().cwd()
    pdb_dir = Path(pdb_dir).resolve()
    os.chdir(pdb_dir)
    command = f'{chimera_executable} --nogui "{ALIGN_SCRIPT_FILE}"'
    print(command)
    subprocess.run(command, shell=True, check=True)
    os.chdir(cwd)


def align_and_export_pdbs_chimerax(
    pdb_dir: str | Path,
    chimera_executable: str = config.CHIMERAX_EXECUTABLE,
):
    """
    Changes the working directory to the pdb_dir, runs the script in there
    """
    cwd = Path().cwd()
    pdb_dir = Path(pdb_dir).resolve()
    os.chdir(pdb_dir)
    command = f'{chimera_executable} --nogui "{ALIGN_EXPORT_SCRIPT_FILE}"'
    print(command)
    subprocess.run(command, shell=True, check=True)
    os.chdir(cwd)


# def align_pdbs_chimera_copy_method(
#     pdb_dir: str | Path,
#     chimera_executable: str = env.CHIMERAX_EXECUTABLE,
#     cleanup: bool = True,
# ):
#     """creates a temporary script file, copies it to the pdb_dir, and runs it.
#     pdb_dir: the directory containing the pdbs to be aligned
#     chimera_executable: the path to the chimera executable
#     """
#     # cwd = Path().cwd()
#     # os.chdir(pdb_dir)
#     pdb_dir = Path(pdb_dir).resolve()
#     shutil.copy(ALIGN_SCRIPT_FILE, pdb_dir)
#     script_path = pdb_dir / ALIGN_SCRIPT_FILE.name
#     command = f'{chimera_executable} --nogui "{script_path}"'
#     print(command)
#     subprocess.run(command, shell=True, check=True)
#     if cleanup:
#         os.remove(script_path)

