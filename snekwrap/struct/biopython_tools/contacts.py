from Bio.PDB import PDBParser, NeighborSearch, Superimposer, Select, Structure  # type: ignore
from Bio.Data import IUPACData
from pathlib import Path
import json
import re
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning

warnings.filterwarnings("ignore", category=PDBConstructionWarning)


def _contacts_to_strings(contacts):
    '''
    This function borrowed/adapted from original FragFold https://github.com/swanss/FragFold
    '''
    return [_get_contact_name(rA, rB) for rA, rB in contacts]


def _get_residue_name(res):
    '''
    This function borrowed/adapted from original FragFold https://github.com/swanss/FragFold
    '''
    return f"{res.get_parent().id}_{res.id[1]}"


def _get_contact_name(resi, resj):
    '''
    This function borrowed/adapted from original FragFold https://github.com/swanss/FragFold
    '''
    # not sure why the order of the residues is important here.
    if resi.get_parent().id <= resj.get_parent().id:
        return f"{_get_residue_name(resi)}-{_get_residue_name(resj)}"
    else:
        return f"{_get_residue_name(resj)}-{_get_residue_name(resi)}"


def _is_interchain_contact(
    res_1,
    res_2,
    chain_group_a: None | set | list = None,
    chain_group_b: None | set | list = None,
):
    """
    returns True if the residues are in different chains or in different chain groups (if defined)
    This function adapted from original FragFold https://github.com/swanss/FragFold
    """
    res1_chain = res_1.get_parent().id
    res2_chain = res_2.get_parent().id
    if chain_group_a is None and chain_group_b is None:
        return res1_chain != res2_chain
    assert (
        chain_group_a is not None and chain_group_b is not None
    ), f"if 1 chain group is defined, you must define the other. chain groups: {chain_group_a=}, {chain_group_b=}"
    if res1_chain in chain_group_a and res2_chain in chain_group_b:
        return True
    if res1_chain in chain_group_b and res2_chain in chain_group_a:
        return True
    return False


def get_interchain_contacts_from_pdb(
    pdb_file: str | Path,
    distance_cutoff: float | int = 4.0,
    chain_group_a: list[str] | None = None,
    chain_group_b: list[str] | None = None,
):
    """
    This function adapted from original FragFold https://github.com/swanss/FragFold
    """
    pdb_file = Path(pdb_file)
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("s", pdb_file)
    contacts = _contacts_to_strings(
        get_interchain_contacts(
            s,
            contact_distance=distance_cutoff,
            chain_group_a=chain_group_a,
            chain_group_b=chain_group_b,
        )
    )
    return contacts


def get_interchain_contacts(
    structure: Structure.Structure,
    contact_distance=4.0,
    chain_group_a=None,
    chain_group_b=None,
):
    """
    This function adapted from original FragFold https://github.com/swanss/FragFold
    """
    ns = NeighborSearch([x for x in structure.get_atoms()])
    nearby_res = ns.search_all(contact_distance, "R")
    contacts = [
        (x, y)
        for x, y in nearby_res  # type: ignore
        if _is_interchain_contact(x, y, chain_group_a, chain_group_b)
    ]
    return contacts


class Residue:

    def __init__(self, res):
        self.pos = res.id[1]
        self.chain = res.get_parent().id
        self.residue = IUPACData.protein_letters_3to1[res.get_resname().capitalize()]

    def __repr__(self):
        return f"{self.chain}.{self.pos} {self.residue}"


class Contact:

    def __init__(self, res1, res2):
        self.res1 = Residue(res1)
        self.res2 = Residue(res2)
        self.bychain = {}
        self.bychain[res1.get_parent().id] = self.res1
        self.bychain[res2.get_parent().id] = self.res2
        # self.res1, self.res2 = sorted([self.res1, self.res2], key=lambda x: (x.chain, x.pos))

    def __repr__(self):
        return f"{self.res1} - {self.res2}"


class StructureContacts:

    def __init__(
        self,
        pdb_path: Path | str,
        contact_distance: float = 5,
        chain_group_a: list[str] | None = None,
        chain_group_b: list[str] | None = None,
    ):
        self.pdb_path = pdb_path
        self.contact_distance = contact_distance
        self.parser = PDBParser(QUIET=True)
        self.chain_group_a = chain_group_a
        self.chain_group_b = chain_group_b
        self.structure = self.parser.get_structure("s", pdb_path)
        self.chains = [chain.id for chain in self.structure.get_chains()]
        self.chain_seqs = self._structure2sequences()
        self.contacts = self._get_contacts()
        self.contact_positions = self._get_contact_positions_per_chain()

    def _structure2sequences(self):
        sequences = {}
        for model in self.structure:
            for chain in model:
                sequences[chain.id] = [[], ""]
                for res in chain:
                    sequences[chain.id][0].append(res.get_id()[1])
                    resname = res.get_resname().capitalize()
                    if resname in IUPACData.protein_letters_3to1:
                        one_letter = IUPACData.protein_letters_3to1[resname]
                    else:
                        raise ValueError(f"Unknown residue name: {resname}")
                    sequences[chain.id][1] += one_letter
        return sequences

    def get_sequence_with_gaps(self, chain_id):
        '''
        Get the sequence of a chain including gaps based on residue positions. 
        gaps should be where residues are missing from the structure.
        '''
        seq_list, seq_str = self.chain_seqs[chain_id]
        if len(seq_list) != len(seq_str):
            raise ValueError(f"Sequence length mismatch for chain {chain_id}")
        seq_with_gaps = []
        last_pos = 0
        for pos, aa in zip(seq_list, seq_str):
            gap_length = pos - last_pos - 1
            if gap_length > 0:
                seq_with_gaps.append("-" * gap_length)
            seq_with_gaps.append(aa)
            last_pos = pos
        return "".join(seq_with_gaps)

    def _get_contact_positions_per_chain(self):
        contacts_per_chain = {}
        for chain_id in self.chains:
            contacts_per_chain[chain_id] = sorted(
                list(set([c.bychain[chain_id].pos for c in self.get_contacts_by_chain(chain_id)]))
            )
        return contacts_per_chain

    def _get_contacts(self):
        """Get contacts between residues in the structure."""
        _contacts = get_interchain_contacts(
            self.structure,
            contact_distance=self.contact_distance,
            chain_group_a=self.chain_group_a,
            chain_group_b=self.chain_group_b,
        )
        return [Contact(i[0], i[1]) for i in _contacts]

    def get_contacts_by_chain(self, chain):
        """Get contacts for a specific chain."""
        contacts_by_chain = []
        for contact in self.contacts:
            if chain in contact.bychain:
                contacts_by_chain.append(contact)
        return contacts_by_chain

    def get_contacts_by_position(self, chain, position):
        """Get contacts for a specific chain and position."""
        contacts_by_position = []
        for contact in self.contacts:
            if contact.res1.pos == position and contact.res1.chain == chain:
                contacts_by_position.append(contact)
            elif contact.res2.pos == position and contact.res2.chain == chain:
                contacts_by_position.append(contact)
        return contacts_by_position
