
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import Align, AlignIO, Seq, SeqIO

from Bio.SeqRecord import SeqRecord
import snekwrap.config as config
import snekwrap.backend.sequence_utils as tools
from typing import Literal


def mafft_align_wrapper(
    input_seqrecord_list: list[SeqRecord],
    mafft_executable: str = config.MAFFT,
    extra_args: str = "",
    n_align_threads: int = 8,
    output_format: str = "dict",
) -> tuple[str, dict | list]:
    # example extra_args: "--retree 1"
    # create temporary file
    temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
    # write seqrecords to temporary file
    SeqIO.write(input_seqrecord_list, temp_file, "fasta")
    temp_file.close()
    # run mafft
    alignment_filename = f"{temp_file.name}-mafft.fa"
    # raise an error if the alignment file already exists. (it won't but just in case)
    if os.path.exists(alignment_filename):
        raise FileExistsError(f"{alignment_filename} already exists")
    else:
        mafft_command = f'{mafft_executable} --thread {n_align_threads} --quiet --anysymbol {extra_args} "{temp_file.name}" > "{alignment_filename}"'
    # print(mafft_command)
    subprocess.run(mafft_command, shell=True, check=True)
    mafft_output = tools.import_fasta(alignment_filename, output_format=output_format)
    # delete temporary file
    os.remove(alignment_filename)
    os.remove(temp_file.name)
    return mafft_command, mafft_output  # type: ignore


def clustal_align_wrapper(
    input_seqrecord_list,
    clustal_executable: str = config.CLUSTALO,
    extra_args: str = "",  # try "--full"
    n_align_threads: int = 8,
    output_type="list",
):
    assert output_type in [
        "list",
        "dict",
        "alignment",
    ], f'`output_type` must be one of ["list", "dict", "alignment"]'
    temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
    SeqIO.write(input_seqrecord_list, temp_file, "fasta")
    temp_file.close()
    alignment_filename = f"{temp_file.name}-clustal.fa"
    if os.path.exists(alignment_filename):
        raise FileExistsError(f"{alignment_filename} already exists")
    clustal_command = f'{clustal_executable} -i "{temp_file.name}" -o "{alignment_filename}" -v --outfmt=fa --threads={n_align_threads} {extra_args}'
    subprocess.run(clustal_command, shell=True, check=True)
    if output_type == "list":
        clustal_output = tools.import_fasta(alignment_filename, output_format="list")
    elif output_type == "dict":
        clustal_output = tools.import_fasta(alignment_filename, output_format="dict")
    else:
        clustal_output = AlignIO.read(alignment_filename, "fasta")
    os.remove(alignment_filename)
    os.remove(temp_file.name)
    return clustal_output


def muscle_align_wrapper(
    input_seqrecord_list: list[SeqRecord],
    muscle_binary: str = config._EXECUTABLES["muscle"],
    output_type: str = "list",
    n_align_threads: int = 8,
) -> list[SeqRecord] | dict[str, SeqRecord] | Align.MultipleSeqAlignment:
    assert output_type in [
        "list",
        "dict",
        "alignment",
    ], f'`output_type` must be one of ["list", "dict", "alignment"]'

    temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
    SeqIO.write(input_seqrecord_list, temp_file, "fasta")
    temp_file.close()
    alignment_filename = f"{temp_file.name}-muscle.fa"
    # raise an error if the alignment file already exists. (it won't but just in case)
    if os.path.exists(alignment_filename):
        raise FileExistsError(f"{alignment_filename} already exists")

    muscle_command = f'{muscle_binary} -super5 "{temp_file.name}" -output "{alignment_filename}" -threads {n_align_threads}'
    subprocess.run(muscle_command, shell=True, check=True)

    if output_type == "list":
        muscle_output = tools.import_fasta(alignment_filename, output_format="list")
    elif output_type == "dict":
        muscle_output = tools.import_fasta(alignment_filename, output_format="dict")
    # elif output_type == "alignment":
    else:
        muscle_output = AlignIO.read(alignment_filename, "fasta")

    # delete temporary file
    os.remove(alignment_filename)
    os.remove(temp_file.name)
    return muscle_output
