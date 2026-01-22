from IPython.display import display
import py3Dmol
from io import StringIO
from Bio.PDB import PDBIO, MMCIFParser, PPBuilder, Structure  # type: ignore
import nglview as nv


def _build_residue_selection_by_name(
    structure: Structure.Structure, residue_names: tuple[str, ...] = ("CYS", "CYX")
) -> str:
    """
    Build an nglview selection string targeting residues by name, grouped per chain.
    """
    selections: list[str] = []
    for chain in structure.get_chains():
        for residue in chain.get_residues():
            if residue.get_resname().strip() not in residue_names:
                continue
            hetfield, resseq, icode = residue.get_id()
            if hetfield.strip() and hetfield != " ":
                continue
            insertion_code = icode.strip()
            resi = f"{resseq}{insertion_code}" if insertion_code else f"{resseq}"
            selections.append(f"(chain {chain.id} and resi {resi})")
    return " or ".join(selections)

def show_structure_with_py3Dmol(structure: Structure.Structure, style: str="cartoon"):
    """
    Displays a Biopython Structure object using py3Dmol in a Jupyter environment.
    """
    # Write structure to a PDB string
    pdb_buf = StringIO()
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_buf)
    pdb_str = pdb_buf.getvalue()

    # Visualize with py3Dmol
    view = py3Dmol.view(width=600, height=400)
    view.addModel(pdb_str, "pdb")
    view.setStyle({style: {}})
    view.zoomTo()
    display(view)


def show_structure_with_nglview(structure: Structure.Structure, style: str="cartoon", highlight_cysteines: bool=False):
    """
    Displays a Biopython Structure object using nglview in a Jupyter environment.
    """
    import nglview as nv

    # Write structure to a PDB string
    pdb_buf = StringIO()
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdb_buf)
    pdb_str = pdb_buf.getvalue()

    # Visualize with nglview
    view = nv.NGLWidget()
    view.add_component(nv.TextStructure(pdb_str), ext="pdb")
    view.clear_representations()
    view.add_representation(style)
    if highlight_cysteines:
        view.add_representation(
            "ball+stick",
            selection="protein AND CYS",
            color="yellow",
        )
    view.center()
    display(view)



        # cys_selection = _build_residue_selection_by_name(structure)
        # if cys_selection:
        #     view.add_representation(
        #         "ball+stick",
        #         selection=cys_selection,
        #         color="yellow",
        #     )
