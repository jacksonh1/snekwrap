from pathlib import Path

from dotenv import load_dotenv
from loguru import logger
from importlib_resources import files
import yaml
import pandas as pd
from dataclasses import dataclass
import os
import copy

_three_to_one_dict = {
    'Ala': 'A',
    'Cys': 'C',
    'Asp': 'D',
    'Glu': 'E',
    'Phe': 'F',
    'Gly': 'G',
    'His': 'H',
    'Ile': 'I',
    'Lys': 'K',
    'Leu': 'L',
    'Met': 'M',
    'Asn': 'N',
    'Pro': 'P',
    'Gln': 'Q',
    'Arg': 'R',
    'Ser': 'S',
    'Thr': 'T',
    'Val': 'V',
    'Trp': 'W',
    'Tyr': 'Y'
}
THREE_TO_ONE_DICT = {k.upper(): v for k, v in _three_to_one_dict.items()}
# Load environment variables from .env file if it exists
load_dotenv(override=True)

# ==============================================================================
# // main project paths
# ==============================================================================
# PROJ_ROOT = Path(__file__).resolve().parents[1]
# logger.info(f"PROJ_ROOT path is: {PROJ_ROOT}")

STANDARD_AMINO_ACIDS = sorted(
    [
        "V",
        "E",
        "K",
        "I",
        "H",
        "L",
        "G",
        "T",
        "M",
        "N",
        "S",
        "P",
        "A",
        "F",
        "W",
        "Y",
        "Q",
        "R",
        "C",
        "D",
    ]
)


# ==============================================================================
# // import executables.yaml file to get the paths of external executables
# ==============================================================================

_default_executables = {
    "colabfold_batch": "colabfold_batch",
    "colabfold_data": "/path/to/colabfold_data",
    "chimerax": "chimerax",
    "mafft": "mafft",
    "clustalo": "clustalo",
    "cd_hit": "cd-hit",
    "USalign": "USalign",
    "muscle": "muscle",
    "proteinmpnn_repo": "../ProteinMPNN/",
    "pottsmpnn_repo": "../PottsMPNN/",
    "rfdiffusion_singularity": "rfd.sif",
    "rfdiffusion_model_dir": "/path/to/rfdiffusion_models",
    # BioEmu related defaults
    # If you want to run BioEmu inside a dedicated conda env, set the env name here
    # or via environment variable BIOEMU_CONDA_ENV. If you have a specific executable
    # path, set bioemu_executable or env var BIOEMU_EXECUTABLE.
    "bioemu_executable": "bioemu",
    "bioemu_conda_env": None,
}
executables_file = [
    i for i in files("snekwrap").iterdir() if i.name == "executables.yaml"
][0]
if executables_file is None:
    print("executables.yaml file not found, assuming all executables are in PATH")
    print("colabfold_data path needs to be set manually.")
    print("Use export COLABFOLD_DATA=/path/to/colabfold_data")
    print("Or set COLABFOLD_DATA in a .env file")
    _EXECUTABLES = copy.deepcopy(_default_executables)
else:
    with executables_file.open() as f:
        _EXECUTABLES = yaml.safe_load(f)
for exe in _default_executables.keys():
    if exe not in _EXECUTABLES:
        _EXECUTABLES[exe] = _default_executables[exe]
COLABFOLD_BATCH = os.environ.get("COLABFOLD_BATCH", _EXECUTABLES["colabfold_batch"])
COLABFOLD_DATA = os.environ.get("COLABFOLD_DATA", _EXECUTABLES["colabfold_data"])
# CHIMERAX = os.environ.get("CHIMERAX", EXECUTABLES["chimerax"])
MAFFT = os.environ.get("MAFFT", _EXECUTABLES["mafft"])
CD_HIT = os.environ.get("CD_HIT", _EXECUTABLES["cd_hit"])
USALIGN = os.environ.get("USALIGN", _EXECUTABLES["USalign"])
MUSCLE = os.environ.get("MUSCLE", _EXECUTABLES["muscle"])
CLUSTALO = os.environ.get("CLUSTALO", _EXECUTABLES["clustalo"])
PROTEINMPNN_REPO = os.environ.get("PROTEINMPNN_REPO", _EXECUTABLES["proteinmpnn_repo"])
POTTSMPNN_REPO = os.environ.get("POTTSMPNN_REPO", _EXECUTABLES["pottsmpnn_repo"])
RFDIFFUSION_SINGULARITY = os.environ.get("RFDIFFUSION_SINGULARITY", _EXECUTABLES["rfdiffusion_singularity"])
RFDIFFUSION_MODEL_DIR = os.environ.get("RFDIFFUSION_MODEL_DIR", _EXECUTABLES["rfdiffusion_model_dir"])

# BioEmu
BIOEMU_EXECUTABLE = os.environ.get("BIOEMU_EXECUTABLE", _EXECUTABLES.get("bioemu_executable", "bioemu"))
# If BIOEMU_CONDA_ENV is unset or empty, the wrapper will run without conda run
BIOEMU_CONDA_ENV = os.environ.get("BIOEMU_CONDA_ENV", _EXECUTABLES.get("bioemu_conda_env"))


COLABFOLD_PDB_PREDICTION_FILENAME_REGEX = r"(?P<name>.+)_\w+_rank_(?P<rank>\d+)_(?P<weights>.+)_model_._seed_\d\d\d.*\.pdb"





