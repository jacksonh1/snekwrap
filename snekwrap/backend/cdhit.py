from pathlib import Path

import snekwrap.backend.sequence_utils as tools
import snekwrap.config as config
import os
import subprocess
import tempfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from typing import Literal



def cd_hit_clstr_parser(clstr_filepath: str | Path) -> dict[str, dict[str, list[str]]]:
    """
    returns dictionary of cluster names and their members
    the dictionary will have two keys for each cluster:
    'all_members' (containing the seq ids for every sequence in the cluster including the representative sequence) and 'representative_seq' (just the id of the representative sequence)
    TODO: might want to not split the sequence id from the ... at the end (removes %ids from cd-hit)
    `clstr_filepath` is the path to the cd-hit output file
    """
    clusters_dict = {}
    cluster = None
    with open(clstr_filepath, "r") as f:
        for line in f:
            if line.startswith(">"):
                cluster = line.split(">")[1].strip()
                clusters_dict[cluster] = {}
                continue
            else:
                if cluster is None:
                    raise ValueError(
                        "Cluster name not found for some of the lines. Check the file format."
                    )
                clusters_dict[cluster].setdefault("all_members", []).append(
                    line.split(">")[1].strip().split("...")[0]
                )
                if line.strip().endswith("*"):
                    clusters_dict[cluster]["representative_seq"] = (
                        line.split(">")[1].strip().split("...")[0]
                    )
    return clusters_dict


def cd_hit_clstr_redefine_cluster_representative_by_keywords(clusters_dict, keywords):
    """
    Redefine cluster representative by keywords.
    :param clusters_dict: dict of clusters
    :param keywords: list of keywords
    :return: dict of clusters
    """
    for cluster in clusters_dict:
        n_found = 0
        found_list = []
        for i in clusters_dict[cluster]["all_members"]:
            if any([keyword in i for keyword in keywords]):
                print(f"found keyword in {cluster}")
                clusters_dict[cluster]["representative_seq"] = i
                n_found += 1
                found_list.append(i)
        assert (
            n_found <= 1
        ), f"found more than one id with {keywords} in {cluster}\n{clusters_dict[cluster]['all_members']}"
    return clusters_dict


def cdhit_clstr_retrieve_representative_sequences(clstr_dict, seqrecord_dict):
    """
    pull out representative seqs defined in cdhit clstr_dict from full seqrecord_dict
    """
    clustered_seq_dict = {}
    for cluster_id in clstr_dict.keys():
        id_i = clstr_dict[cluster_id]["representative_seq"]
        # id_i = re.findall(r'\d+\_\d\:.+$', rep_i)[0]
        clustered_seq_dict[id_i] = seqrecord_dict[id_i]
    return clustered_seq_dict


def cdhit_minidriver(seqrecords_2_cluster_list, repr_id_keywords):
    """
    in: list of seqrecords
    out: dict of clustered seqrecords
    for getting LDOs first and then clustering
    """
    seqrecords_2_cluster_dict = {
        seqrecord.id: seqrecord for seqrecord in seqrecords_2_cluster_list
    }
    _, _, cdhit_clstr_dict = cd_hit_wrapper(seqrecords_2_cluster_list)
    cdhit_clstr_dict = cd_hit_clstr_redefine_cluster_representative_by_keywords(
        cdhit_clstr_dict, keywords=repr_id_keywords
    )
    sequences_clustered_OG_dict = cdhit_clstr_retrieve_representative_sequences(
        cdhit_clstr_dict, seqrecords_2_cluster_dict
    )
    return sequences_clustered_OG_dict


def cd_hit_wrapper(
    input_seqrecord_list: list[SeqRecord],
    cd_hit_executable: str = config.CD_HIT,
    extra_args: str = "",
) -> tuple[str, dict[str, SeqRecord], dict[str, dict[str, list[str]]]]:

    # create temporary file
    temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
    # write seqrecords to temporary file
    SeqIO.write(input_seqrecord_list, temp_file, "fasta")
    temp_file.close()
    clustered_seqs_filename = f"{temp_file.name}-cdhit.fa"
    # raise an error if the alignment file already exists. (it won't but just in case)
    if os.path.exists(clustered_seqs_filename):
        raise FileExistsError(f"{clustered_seqs_filename} already exists")

    clustered_seqs_clusters_filename = clustered_seqs_filename + ".clstr"
    command = f"{cd_hit_executable} -i {temp_file.name} -o {clustered_seqs_filename} -M 0 -d 0 -g 1 {extra_args}"
    subprocess.run(command, shell=True, check=True)

    output_clstrs_dict = cd_hit_clstr_parser(
        clustered_seqs_clusters_filename
    )

    output = tools.import_fasta(clustered_seqs_filename, output_format="dict")
    # delete temporary file
    os.remove(clustered_seqs_filename)
    os.remove(clustered_seqs_clusters_filename)
    os.remove(temp_file.name)
    return command, output, output_clstrs_dict  # type: ignore