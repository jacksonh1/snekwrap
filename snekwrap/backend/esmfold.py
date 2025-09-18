from transformers import AutoTokenizer, EsmForProteinFolding
import torch


def predict_structure(sequence, output_file="structure.pdb"):
    """
    Predict the structure of a protein sequence using ESMFold and save it to a PDB file.

    Args:
        sequence (str): The amino acid sequence of the protein.
        output_file (str): The path to the output PDB file.
    """
    assert isinstance(sequence, str), "Sequence must be a string"
    assert all(
        aa in "ACDEFGHIKLMNPQRSTVWY" for aa in sequence
    ), "Sequence contains invalid amino acids"
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=False)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)
    model.esm = model.esm.float()
    # Configure chunking for memory efficiency
    # model.trunk.set_chunk_size(64)
    with torch.no_grad():
        output = model.infer_pdb(sequence)

    # Save result
    with open(output_file, "w") as f:
        f.write(output)