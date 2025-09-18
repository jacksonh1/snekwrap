import re
import json
from pathlib import Path
from snekwrap import config
from typing import Literal
import subprocess
import os
import snekwrap.backend.sequence_utils as seq_utils



def colabfold_batch_wrapper(
    input_file_or_directory: str | Path,
    output_dir: str | Path,
    weights: Literal[
        "alphafold2",
        "alphafold2_ptm",
        "alphafold2_multimer_v1",
        "alphafold2_multimer_v2",
        "alphafold2_multimer_v3",
        "deepfold_v1",
    ] = "alphafold2_ptm",
    pairmode: Literal[
        "unpaired",
        "paired",
        "unpaired_paired",
    ] = "unpaired",
    colabfold_executable: str | Path = config.COLABFOLD_BATCH,
    colabfold_data: str | Path = config.COLABFOLD_DATA,
    extra_args: str | Path = "",
    run = True,
    # gpu_number: int = 0,
) -> str:
    # subprocess.run(
    #     "AVAILABLE_GPU=$(nvidia-smi --query-gpu=index,memory.used --format=csv,noheader | sort -n -k2 | head -n1 | cut -d, -f1)",
    #     shell=True,
    #     check=True,
    # )
    # subprocess.run("export CUDA_VISIBLE_DEVICES=$AVAILABLE_GPU", shell=True, check=True)
    # subprocess.run(f"export CUDA_VISIBLE_DEVICES={gpu_number}", shell=True, check=True)
    # colab_command = f"{colabfold_executable} --model-type alphafold2_ptm --amber --num-relax 5 --use-gpu-relax --data '{colabfold_data}' --pair-mode unpaired {input_file_or_directory} {output_dir}"
    # colab_command = f"export CUDA_VISIBLE_DEVICES={gpu_number}; {colabfold_executable} --model-type {weights} --data '{colabfold_data}' --pair-mode {pairmode} {extra_args} {input_file_or_directory} {output_dir}"
    colab_command = f"{colabfold_executable} --model-type {weights} --data '{colabfold_data}' --pair-mode {pairmode} {extra_args} '{input_file_or_directory}' '{output_dir}'"
    if run:
        subprocess.run(colab_command, shell=True, check=True)
        return colab_command
    else:
        return colab_command


def colabfold_batch_MSA_wrapper(
    input_file: str | Path,
    output_dir: str | Path,
    colabfold_executable: str = config.COLABFOLD_BATCH,
    colabfold_data: str = config.COLABFOLD_DATA,
):
    # subprocess.run("export MPLBACKEND=Agg", shell=True, check=True)
    os.environ["JAX_PLATFORMS"] = "cpu"
    colab_command = f'{colabfold_executable} --msa-only --data "{colabfold_data}" --msa-only "{input_file}" "{output_dir}"'
    subprocess.run(colab_command, shell=True, check=True)



def get_colabfold_msa(
    fasta_file: str | Path, msa_cache_dir: str | Path, **kwargs
) -> Path:
    fasta_file = Path(fasta_file)
    seqs = seq_utils.import_fasta(fasta_file)
    if len(seqs) != 1:
        raise ValueError(
            f"Input fasta file {fasta_file} must contain exactly one sequence."
        )
    seqid = str(seqs[0].id)  # type: ignore
    # filestem = fasta_file.stem
    msa_file = Path(msa_cache_dir) / f"{seqid}.a3m"
    if not msa_file.exists():
        # Generate MSA using colabfold or any other method
        # This is a placeholder for the actual MSA generation logic
        colabfold_batch_MSA_wrapper(
            input_file=fasta_file,
            output_dir=msa_cache_dir,
            **kwargs,  # pass any additional arguments needed for MSA generation
        )
        if not msa_file.exists():
            raise FileNotFoundError(
                f"MSA file {msa_file} was not found despite attempting download. \n check {msa_cache_dir} and {fasta_file} sequence id"
            )
    else:
        print(
            f"MSA file ({msa_file}) with same name found in {msa_cache_dir}, skipping download."
        )
    return msa_file




def get_colabfold_scores(pdb_file: str | Path):
    pdb_file = Path(pdb_file)
    score_file = pdb_file.stem.replace("_unrelaxed", "_scores") + ".json"
    score_file = pdb_file.parent / score_file
    with open(score_file) as f:
        scores = json.load(f)
    return score_file, scores


def colabfold_pdb_filename_2_score_filename(pdb_filename):
    return pdb_filename.stem.replace("_unrelaxed", "_scores") + ".json"


def parse_rank_from_pdb_filename(
    pdb_filename, pattern=r".+_unrelaxed_rank_(\d\d\d)_.+"
):
    """
    extracts the rank from the name of the `pdb_file` using the provided regex pattern.

    Parameters
    ----------
    pdb_filename : str
        name of the pdb file
    pattern : regexp, optional
        regex pattern to extract structure rank from the `pdb_filename`,
        by default r".+_unrelaxed_rank_(\\d\\d\\d)_.+"
    """
    filename_pattern = re.compile(pattern)
    match = filename_pattern.match(pdb_filename)
    assert (
        match is not None
    ), f"File name {pdb_filename} does not match pattern {pattern}"
    rank = match.group(1)
    return int(rank)


def parse_pdb_filename_general(pdb_filename, pattern=config.COLABFOLD_PDB_PREDICTION_FILENAME_REGEX):
    """
    extracts the rank, filename and weights from the name of the `pdb_file` using the provided regex pattern.

    Parameters
    ----------
    pdb_filename : str
        name of the pdb file
    pattern : regexp, optional
        regex pattern to extract structure rank from the `pdb_filename`,
        by default r"(?P<name>.+)_\w+_rank_(?P<rank>\d+)_(?P<weights>.+)_model_._seed_\d\d\d.*\.pdb"
    """
    pdb_filename = Path(pdb_filename)
    p = re.compile(pattern)
    m = p.match(pdb_filename.name)
    if m is None:
        raise ValueError(f"Filename {pdb_filename.name} does not match pattern {pattern}")
    return m.groupdict()
