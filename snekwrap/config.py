from pathlib import Path

from dotenv import load_dotenv
from loguru import logger
from importlib_resources import files
import yaml
import pandas as pd
from dataclasses import dataclass
import os
import copy

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
    "cd_hit": "cd-hit",
    "USalign": "USalign",
    "muscle": "muscle",
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
CLUSTALO = os.environ.get("CLUSTALO", _EXECUTABLES.get("clustalo", "clustalo"))




COLABFOLD_PDB_PREDICTION_FILENAME_REGEX = r"(?P<name>.+)_\w+_rank_(?P<rank>\d+)_(?P<weights>.+)_model_._seed_\d\d\d.*\.pdb"





