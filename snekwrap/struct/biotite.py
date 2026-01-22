import biotite.structure.io.pdb as pdb
import biotite.structure as struc
import structutils.config as config
import math

def validate_sequence_numbering(sequence: str, atom_array):
    """
    Validate that each present residue in the atom array matches the residue
    at the same 1-based position in the provided sequence.

    Parameters
    - sequence: str, full sequence with positions starting at 1
    - atom_array: biotite.structure.AtomArray for a single chain

    Returns
    - List of dicts describing mismatches or out-of-range indices.
      Empty list means all present residues match.
    """
    seq = sequence.strip().upper()
    mismatches = []
    seen_res_ids = set()
    for res_id, res_name in zip(atom_array.res_id, atom_array.res_name):
        # Only check each residue number once
        if res_id in seen_res_ids:
            continue
        seen_res_ids.add(int(res_id))

        aa1 = config.THREE_TO_ONE_DICT.get(res_name.upper())
        # Skip non-standard residues not in the mapping
        if aa1 is None:
            continue
        idx = int(res_id) - 1  # convert to 0-based index
        if idx < 0 or idx >= len(seq):
            print(f"warning: residue ID {res_id} ({res_name}) is out of range for sequence of length {len(seq)}")
            continue
        if seq[idx] != aa1:
            print(f"Residue mismatch at position {res_id}: structure has {aa1}, sequence has {seq[idx]}")
            mismatches.append({
                "res_id": res_id,
                "structure_residue": aa1,
                "sequence_residue": seq[idx]
            })
    print(f"Total mismatches found: {len(mismatches)}")
    return mismatches


def get_sequence_with_gaps(atom_array):
    """
    Generate a sequence string with gaps ('-') for missing residues
    based on the residue IDs in the atom array.

    Parameters
    - atom_array: biotite.structure.AtomArray for a single chain

    Returns
    - str: sequence with gaps for missing residues
    """
    chain_ = set(atom_array.chain_id)
    if len(chain_) > 1:
        raise ValueError("Atom array must contain only one chain.")
    ca_atoms = atom_array[atom_array.atom_name == "CA"]
    res_ids = sorted(set(int(rid) for rid in ca_atoms.res_id))
    min_res_id = min(res_ids)
    max_res_id = max(res_ids)
    seq_with_gaps = []
    res_names = set(ca_atoms.res_name)
    for res_name in res_names:
        if res_name.upper() not in config.THREE_TO_ONE_DICT:
            print(f"warning: non-standard residue {res_name} found in chain {chain_}. It will be treated as gap.")
    res_id_to_aa = {int(rid): config.THREE_TO_ONE_DICT.get(res_name.upper())
                    for rid, res_name in zip(ca_atoms.res_id, ca_atoms.res_name)}
    for res_id in range(min_res_id, max_res_id + 1):
        aa = res_id_to_aa.get(res_id, '-')
        seq_with_gaps.append(aa)
    return ''.join(seq_with_gaps)


def get_radius_of_gyration(atom_array):
    """
    Calculate the radius of gyration for the given atom array.

    Parameters
    - atom_array: biotite.structure.AtomArray

    Returns
    - float: radius of gyration in angstroms
    """
    coords = atom_array.coord
    center_of_mass = coords.mean(axis=0)
    squared_distances = ((coords - center_of_mass) ** 2).sum(axis=1)
    rg = math.sqrt(squared_distances.mean())
    return rg