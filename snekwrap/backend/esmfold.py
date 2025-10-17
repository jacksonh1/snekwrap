from transformers import AutoTokenizer, EsmForProteinFolding
import torch
from pathlib import Path


def predict_structure(sequence, output_file: str | Path = "structure.pdb"):
    """
    Predict the structure of a protein sequence using ESMFold and save it to a PDB file.

    Args:
        sequence (str): The amino acid sequence of the protein.
        output_file (str): The path to the output PDB file.
    """
    assert isinstance(sequence, str), "Sequence must be a string"
    if ":" in sequence:
        sequence = sequence.replace(":", "X" * 50)
    # assert all(
    #     aa in "ACDEFGHIKLMNPQRSTVWY" for aa in sequence
    # ), "Sequence contains invalid amino acids"
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=False)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)
    model.esm = model.esm.float()
    # Configure chunking for memory efficiency
    model.trunk.set_chunk_size(64)
    with torch.no_grad():
        output = model.infer_pdb(sequence)

    # Save result
    with open(output_file, "w") as f:
        f.write(output)


from transformers import AutoTokenizer, EsmForProteinFolding
import torch
from pathlib import Path
from typing import List, Union
import re


def predict_multichain_with_linkers(
    sequences: Union[str, List[str]],
    output_file: Union[str, Path] = "complex.pdb",
    linker_length: int = 25,
):
    """
    Predict multichain complex using glycine linkers instead of colons
    """

    # Handle input format
    if isinstance(sequences, list):
        chain_sequences = sequences
    elif isinstance(sequences, str):
        if ":" in sequences:
            chain_sequences = sequences.split(":")
        else:
            chain_sequences = [sequences]
    # Create linker
    linker = "G" * linker_length
    print(f"Using {linker_length}-residue glycine linker: {linker}")

    # Combine sequences with linkers
    if len(chain_sequences) > 1:
        combined_sequence = linker.join(chain_sequences)
        chain_starts = []
        current_pos = 0

        for i, seq in enumerate(chain_sequences):
            chain_starts.append(current_pos)
            current_pos += len(seq)
            if i < len(chain_sequences) - 1:  # Add linker length except for last chain
                current_pos += linker_length

        print(f"Combined sequence length: {len(combined_sequence)} residues")
        print(f"Chain start positions: {chain_starts}")
    else:
        combined_sequence = chain_sequences[0]
        chain_starts = [0]

    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=False)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)
    model.esm = model.esm.float()
    model.trunk.set_chunk_size(64)
    with torch.no_grad():
        output = model.infer_pdb(combined_sequence)
    # Post-process PDB to assign proper chain IDs and remove linkers
    if len(chain_sequences) > 1:
        processed_pdb = process_multichain_pdb(
            output, chain_sequences, chain_starts, linker_length
        )
    else:
        processed_pdb = output

    # Save result
    with open(output_file, "w") as f:
        f.write(processed_pdb)


def process_multichain_pdb(
    pdb_content: str, chain_sequences: List[str], chain_starts: List[int], linker_length: int
) -> str:
    """
    Process PDB to assign chain IDs and remove linker residues
    """
    lines = pdb_content.split("\n")
    processed_lines = []
    chain_letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # Get residue number (1-indexed in PDB)
            res_num = int(line[22:26].strip())
            seq_pos = res_num - 1  # Convert to 0-indexed

            # Determine which chain this residue belongs to
            current_chain = 0
            adjusted_res_num = seq_pos

            for i, start_pos in enumerate(chain_starts):
                chain_end = start_pos + len(chain_sequences[i])
                if i < len(chain_starts) - 1:
                    next_chain_start = chain_starts[i + 1]
                    # Check if in linker region (remove these)
                    if chain_end <= seq_pos < next_chain_start:
                        # Skip linker residues
                        continue

                if start_pos <= seq_pos < chain_end:
                    current_chain = i
                    adjusted_res_num = seq_pos - start_pos + 1
                    break
                elif seq_pos >= chain_end:
                    # Adjust for linkers in previous chains
                    adjusted_res_num -= linker_length * (i + 1)

            # Skip if this is a linker residue
            skip_residue = False
            for i in range(len(chain_sequences) - 1):
                linker_start = chain_starts[i] + len(chain_sequences[i])
                linker_end = chain_starts[i + 1]
                if linker_start <= seq_pos < linker_end:
                    skip_residue = True
                    break

            if skip_residue:
                continue

            # Assign chain ID and adjust residue number
            chain_id = chain_letters[current_chain]
            new_line = line[:21] + chain_id + f"{adjusted_res_num:4d}" + line[26:]
            processed_lines.append(new_line)
        else:
            processed_lines.append(line)

    return "\n".join(processed_lines)
