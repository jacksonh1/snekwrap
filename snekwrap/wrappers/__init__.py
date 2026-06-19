import importlib
from typing import TYPE_CHECKING

from loguru import logger

# Wrappers are imported lazily (PEP 562): each subwrapper loads on first attribute
# access rather than at `import snekwrap` time. Several of these pull in torch /
# transformers, so eager-importing them all made `import snekwrap` slow even when
# only one wrapper was needed.
_SUBMODULES = {
    "cdhit": "snekwrap.wrappers.cdhit",
    "colabfold": "snekwrap.wrappers.colabfold",
    "esmfold": "snekwrap.wrappers.esmfold",
    "MSAs": "snekwrap.wrappers.MSAs",
    "protein_mpnn": "snekwrap.wrappers.proteinmpnn.protein_mpnn",
    "RFDiffusion": "snekwrap.wrappers.rfdiff.RFDiffusion",
    "potts_mpnn": "snekwrap.wrappers.pottsmpnn.potts_mpnn",
}

__all__ = list(_SUBMODULES)


def __getattr__(name: str):
    """Lazily import and cache a wrapper submodule on first access (PEP 562)."""
    if name in _SUBMODULES:
        logger.info(f"importing {name}")
        module = importlib.import_module(_SUBMODULES[name])
        globals()[name] = module  # cache so subsequent access skips this hook
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return __all__


# Resolved statically by type checkers / language servers for autocomplete and
# go-to-definition; never executed at runtime (so it stays lazy).
if TYPE_CHECKING:
    import snekwrap.wrappers.cdhit as cdhit
    import snekwrap.wrappers.colabfold as colabfold
    import snekwrap.wrappers.esmfold as esmfold
    import snekwrap.wrappers.MSAs as MSAs
    import snekwrap.wrappers.proteinmpnn.protein_mpnn as protein_mpnn
    import snekwrap.wrappers.pottsmpnn.potts_mpnn as potts_mpnn
    import snekwrap.wrappers.rfdiff.RFDiffusion as RFDiffusion
