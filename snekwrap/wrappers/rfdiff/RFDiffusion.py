# %%
import re
import json
from pathlib import Path
from snekwrap import config
from typing import Literal
import subprocess
import os
import snekwrap.seq.seqtools as seq_utils
from biopandas.pdb import PandasPdb
import numpy as np
import re
import tempfile
import shutil
import snekwrap.wrappers.rfdiff.colabdesign_utils as colabdesign_utils
from Bio.PDB import PDBParser, MMCIFParser, PPBuilder # type: ignore
from Bio.SeqIO.PdbIO import CifSeqresIterator
from Bio import PDB
from Bio.PDB import PDBIO
from Bio.PDB import Select

CHAIN_OPTIONS = set(list("ABCDEFGHIJKLMNOPQRSTUVWXYZ"))


def pdb_to_chain_position_dict_biopandas(pdb_file):
    ppdb = PandasPdb()
    p = ppdb.read_pdb(pdb_file)
    chain_res_dict = (
        p.df['ATOM'][['chain_id', 'residue_number']]
        .drop_duplicates()
        .groupby('chain_id')['residue_number']
        .apply(list)
        .to_dict()
    )
    return chain_res_dict


def fix_range_4_missing_residues(contig_str, pdb_chain_position_dict):
    pattern = re.compile(r'([A-Za-z]+)(\d+)-(\d+)$')
    if pattern.match(contig_str):
        chain_id, start, end = pattern.match(contig_str).groups()
        start = int(start)
        end = int(end)
    else:
        raise ValueError(f"Invalid contig format for repairing. Expected format: 'X1-20' or 'Y76-94', etc., got: '{contig_str}'")
    if start > end:
        raise ValueError(f"Start position {start} cannot be greater than end position {end} in contig '{contig_str}'")
    if chain_id not in pdb_chain_position_dict:
        raise ValueError(f"Chain ID '{chain_id}' from contig '{contig_str}' not found in PDB data.")
    positions = sorted(pdb_chain_position_dict[chain_id])
    new_ranges = []
    current_start = start
    current_end = start - 1
    if end > positions[-1]:
        end = positions[-1]
    in_range = False
    for pos in positions:
        if pos >= start and pos <= end:
            in_range = True
            if pos == current_end + 1:
                current_end = pos
            elif current_start == current_end + 1:
                current_start = pos
                current_end = pos
            elif current_start < current_end:
                new_ranges.append(f"{chain_id}{current_start}-{current_end}")
                current_start = pos
                current_end = pos
            else:
                raise ValueError("Unexpected condition encountered while processing positions")
            if pos == end:
                new_ranges.append(f"{chain_id}{current_start}-{current_end}")
                in_range = False
        elif in_range:
            new_ranges.append(f"{chain_id}{current_start}-{current_end}")
            in_range = False
    return new_ranges


def get_unique_chain_set_in_contig(contig):
    pattern = re.compile(r'([A-Za-z]+)(\d+)-(\d+)$')
    chain_set = set()
    for part in contig.split("/"):
        if pattern.match(part):
            chain_id, start, end = pattern.match(part).groups()
            chain_set.add(chain_id)
    return chain_set


def get_chains_from_pdb(path):
    parser = PDBParser()
    structure = parser.get_structure(Path(path).stem, path)
    chain_id_list = [str(c.id) for c in structure.get_chains()]
    return chain_id_list


def import_structure_biopython(pdb_file):
    pdb_file = Path(pdb_file)
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    return structure


def remove_chain(structure, chain_id):
    for model in structure:
        chain_to_remove = None
        for chain in model:
            if chain.id == chain_id:
                chain_to_remove = chain
                break
        if chain_to_remove:
            model.detach_child(chain_to_remove.id)
    return structure


def save_structure_to_file(structure, output_path):
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(output_path))
    return output_path


class ChainRm(Select):

    def __init__(self, chainlist):
        self.chainlist = chainlist

    def accept_chain(self, chain):
        if chain.id in self.chainlist:
            return 0
        else:
            return 1



def combine_chains(input_pdb, output_pdb, contig):
    pattern = re.compile(r'([A-Za-z]+)(\d+)-(\d+)$')
    taken_chains = set(get_chains_from_pdb(input_pdb))
    available_chains = CHAIN_OPTIONS - taken_chains
    new_chain_id = list(available_chains)[0]
    chains_to_remove = []
    new_contig_parts = []
    s = import_structure_biopython(input_pdb)[0]
    for part in contig.split("/"):
        if not pattern.match(part):
            new_contig_parts.append(part)
            continue
        chain_id, start, end = pattern.match(part).groups()
        start = int(start)
        end = int(end)
        chains = [str(c.id) for c in s.get_chains()]
        if new_chain_id in chains:
            chain_end = max([r.get_id()[1] for r in s[new_chain_id].get_residues()])
            new_start = chain_end + 1
            new_end = new_start + (end - start)
        else:
            new_start = 1
            new_end = new_start + (end - start)
            # add new chain to structure
            new_chain = PDB.Chain.Chain(new_chain_id)
            s.add(new_chain)
        # move residues from chain_id new_start to new_end to new_chain_id
        for r in s[chain_id].get_residues():
            if r.get_id()[1] >= start and r.get_id()[1] <= end:
                r_copy = r.copy()
                r_copy.id = (' ', new_start + (r.get_id()[1] - start), ' ')
                s[new_chain_id].add(r_copy)
        new_contig_parts.append(f"{new_chain_id}{new_start}-{new_end}")
        chains_to_remove.append(chain_id)
        s = remove_chain(s, chain_id)
    io = PDBIO()
    io.set_structure(s)
    io.save(output_pdb, ChainRm(chains_to_remove))
    return "/".join(new_contig_parts)


def preprocess_pdb_for_fusion(contig_str, pdb_file):
    contigs = contig_str.split(" ")
    new_contigs = []
    pdb_file = Path(pdb_file)
    temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".pdb")
    shutil.copy(pdb_file, temp_file.name)
    temp_file.close()
    chain_change_map = {}
    for c in contigs:
        chain_set = get_unique_chain_set_in_contig(c)
        if len(chain_set) > 1:
            new_c = combine_chains(temp_file.name, temp_file.name, c)
            new_contigs.append(new_c)
            if c in chain_change_map:
                raise ValueError(f"Contig '{c}' has already been processed for chain combination. Something is very wrong.")
            chain_change_map[c] = new_c
        else:
            new_contigs.append(c)
    return temp_file, " ".join(new_contigs), chain_change_map


def convert_chain_change_map_2_sane_dict(chain_change_map):
    pattern = re.compile(r'([A-Za-z]+)(\d+)-(\d+)$')
    mapping = {}
    for k, v in chain_change_map.items():
        key_str = k.strip('/0')
        value_str = v.strip('/0')
        key_parts = key_str.split("/")
        value_parts = value_str.split("/")
        for key_part, value_part in zip(key_parts, value_parts):
            if pattern.match(key_part):
                old_chain_id, old_start, old_end = pattern.match(key_part).groups()
                old_start = int(old_start)
                old_end = int(old_end)
                new_chain_id, new_start, new_end = pattern.match(value_part).groups()
                new_start = int(new_start)
                new_end = int(new_end)
                if old_chain_id in mapping:
                    raise ValueError(f"Chain ID '{old_chain_id}' appears multiple times in chain change map. Cannot process.")
                mapping[old_chain_id] = {}
                mapping[old_chain_id]['new_chain_id'] = new_chain_id
                mapping[old_chain_id]['residue_offset'] = new_start - old_start
                mapping[old_chain_id]['old_range'] = (old_start, old_end)
                mapping[old_chain_id]['new_range'] = (new_start, new_end)
    return mapping
# r"([A-Za-z]+)(\d+)-(\d+)"





def fix_contigs_avoid_missing_residues(contig_str, pdb_file):
    contigs = contig_str.split(" ")
    chain_dict = pdb_to_chain_position_dict_biopandas(pdb_file)
    new_contigs = []
    for c in contigs:
        new_contig = []
        if "-" not in c and re.match(r'^\d+$', c): #all([not i.isalpha() for i in c]):  # re.match(r'^\d+$', c):
            new_contigs.append(f"{c}-{c}")
            continue
        chain_set = get_unique_chain_set_in_contig(c)
        if len(chain_set) > 1:
            # preprocess pdb
            pass
        for part in c.split("/"):
            if re.match(r'^[A-Za-z]+\d+-\d+$', part):
                new_part = "/".join(fix_range_4_missing_residues(part, chain_dict))
            else:
                new_part = part
            new_contig.append(new_part)
        new_contigs.append("/".join(new_contig))
    return " ".join(new_contigs)


def generate_filename_stem_counter(output_dir: Path, original_filename: Path) -> str:
    stem = original_filename.stem
    suffix = original_filename.suffix or ".pdb"
    candidate = original_filename.name
    counter = 1
    while (output_dir / candidate).exists():
        candidate = f"{stem}_{counter}{suffix}"
        counter += 1
    return candidate


def run_rfdiffusion_singularity(
    input_pdb_file: str | Path,
    output_dir: str | Path,
    output_prefix: str = "rfdiff_design",
    singularity_exec_file: str | Path = config.RFDIFFUSION_SINGULARITY,
    model_dir: str | Path = config.RFDIFFUSION_MODEL_DIR,
    num_designs: int = 8,
    contig_str: str = "",
    # diffuser_steps: int = 50,
    hotspot_residues: str = "",
    write_trajectory: bool = False,
    extra_args: str | Path = "",
    fix_contigs: bool = True,
    run: bool = True,
) -> str:
    '''
    # X1-20/2-5/Y76-94/0 B1-107
    '''
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    original_input_pdb_file = Path(input_pdb_file)
    _processed_input_pdb_filename = Path(original_input_pdb_file.stem + '_rfdiff_input.pdb')
    original_contig_str = contig_str
    chain_change_map = {}
    if fix_contigs:
        print("pre-fusion contig:", contig_str)
        temp_file, contig_str, chain_change_map = preprocess_pdb_for_fusion(contig_str, original_input_pdb_file)# temp file is created but not copied to output dir
        print(f"Post-fusion contig: {contig_str}")
        contig_str = fix_contigs_avoid_missing_residues(contig_str, temp_file.name)
        print(f"Post-missing-residue-fix contig: {contig_str}")
        processed_input_pdb_filename = generate_filename_stem_counter(output_dir, _processed_input_pdb_filename)
        input_pdb_file = output_dir / processed_input_pdb_filename
    else:
        temp_file = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".pdb")
        shutil.copy(original_input_pdb_file.resolve(), temp_file.name)
        temp_file.close()
        processed_input_pdb_filename = generate_filename_stem_counter(output_dir, _processed_input_pdb_filename)
        input_pdb_file = output_dir / processed_input_pdb_filename
    input_pdb_dir = input_pdb_file.parent
    schedule_path = output_dir / "schedules"
    schedule_path.mkdir(parents=True, exist_ok=True)
    rfdiffusion_command = f"""singularity run --nv \
    --writable-tmpfs \
    --pwd /app/RFdiffusion \
    --bind {input_pdb_dir}:$HOME/inputs \
    --bind {output_dir}:$HOME/outputs \
    --bind {schedule_path}:$HOME/schedules \
    {singularity_exec_file} \
    inference.output_prefix=$HOME/outputs/{output_prefix} \
    inference.model_directory_path=/app/models \
    inference.schedule_directory_path=$HOME/schedules \
    inference.input_pdb=$HOME/inputs/{processed_input_pdb_filename} \
    inference.num_designs={num_designs} \
    'contigmap.contigs=[{contig_str}]' \
    inference.write_trajectory={str(write_trajectory).lower()} \
    {extra_args}"""
    produced_pdbs = [output_dir / f"{output_prefix}_{i}.pdb" for i in range(num_designs)]
    arguments = {
        "unprocessed_input_pdb_file": str(original_input_pdb_file.resolve()),
        "processing-chain_changes": chain_change_map,
        "processing-chain_change_conversion_map": convert_chain_change_map_2_sane_dict(chain_change_map),
        "input_pdb_file": str(input_pdb_file),
        "output_dir": str(output_dir),
        "output_prefix": output_prefix,
        "output_pdbs": [str(i) for i in produced_pdbs],
        "singularity_exec_file": str(singularity_exec_file),
        # "model_dir": str(model_dir),
        "num_designs": num_designs,
        "original_contig_str": original_contig_str,
        "contig_str": contig_str,
        # "diffuser_steps": diffuser_steps,
        "hotspot_residues": hotspot_residues,
        "write_trajectory": write_trajectory,
        "extra_args": str(extra_args),
        "fix_contigs": fix_contigs,
    }
    # --bind {model_dir}:$HOME/models \
    # inference.model_directory_path=$HOME/models \
    # ppi.hotspot_residues='[E64,E88,E96]' \
    # diffuser.T={diffuser_steps} # DON'T MESS WITH THIS \
    if all(p.exists() for p in produced_pdbs):
        print("All output PDBs already exist, skipping RFDiffusion run.")
        return rfdiffusion_command, arguments
    rfdiff_params_path = output_dir / "rfdiff_params.json"
    if not rfdiff_params_path.exists():
        existing_arguments = []
    else:
        with open(rfdiff_params_path, "r") as f:
            existing_arguments = json.load(f)
    existing_arguments.append(arguments)
    with open(rfdiff_params_path, "w") as f:
        json.dump(existing_arguments, f, indent=4)
    rfdiff_command_path = output_dir / "rfdiff_command.sh"
    if not rfdiff_command_path.exists():
        rfdiff_command_path.touch()
    with open(rfdiff_command_path, "a") as f:
        f.write(rfdiffusion_command + "\n")
    if run:
        shutil.copy(temp_file.name, input_pdb_file) # copy temp file to output dir
        os.remove(temp_file.name)
        subprocess.run(rfdiffusion_command, shell=True, check=True)
        # contigs = contig_str.replace('/0 ', ' ').split(" ")
        # contigs = contig_str.split(" ")
        # for pdb in produced_pdbs:
        #     with open(pdb, "r") as f:
        #         pdb_str = f.read()
        #     fixed_pdb_str = colabdesign_utils.fix_pdb(pdb_str, contigs)
        #     with open(pdb, "w") as f:
        #         f.write(fixed_pdb_str)
        return rfdiffusion_command, arguments
    else:
        os.remove(temp_file.name)
        return rfdiffusion_command, arguments
# def run_rfdiffusion_docker(


def run_rfdiffusion_docker(
    input_pdb_file: str | Path,
    output_dir: str | Path,
    output_prefix: str = "rfdiff_design",
    model_dir: str | Path = config.RFDIFFUSION_MODEL_DIR,
    num_designs: int = 8,
    contig_str: str = "",
    # diffuser_steps: int = 50,
    hotspot_residues: str = "",
    write_trajectory: bool = False,
    extra_args: str | Path = "",
    fix_contigs: bool = True,
    run: bool = True,
) -> str:
    '''
    # X1-20/2-5/Y76-94/0 B1-107
    '''
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if fix_contigs:
        temp_file, contig_str = preprocess_pdb_for_fusion(contig_str, input_pdb_file)
        print(f"Post-fusion contig: {contig_str}")
        contig_str = fix_contigs_avoid_missing_residues(contig_str, temp_file.name)
        print(f"Post-missing-residue-fix contig: {contig_str}")
        shutil.copy(temp_file.name, output_dir / "input.pdb")
        input_pdb_file = output_dir / "input.pdb"
    else:
        input_pdb_file = Path(input_pdb_file).resolve()
    input_pdb_dir = input_pdb_file.parent
    input_pdb_filename = input_pdb_file.name
    print(input_pdb_dir)
    print(input_pdb_filename)
    schedule_path = output_dir / "schedules"
    schedule_path.mkdir(parents=True, exist_ok=True)
    rfdiffusion_command = f"""docker run -it --rm --gpus all \
    -v {model_dir}:$HOME/models \
    -v {input_pdb_dir}:$HOME/inputs \
    -v {output_dir}:$HOME/outputs \
    -v {schedule_path}:$HOME/schedules \
    rfdiffusion \
    inference.output_prefix=$HOME/outputs/{output_prefix} \
    inference.model_directory_path=$HOME/models \
    inference.schedule_directory_path=$HOME/schedules \
    inference.input_pdb=$HOME/inputs/{input_pdb_filename} \
    inference.num_designs={num_designs} \
    'contigmap.contigs=[{contig_str}]' \
    inference.write_trajectory={str(write_trajectory).lower()} \
    {extra_args}"""
    # ppi.hotspot_residues='[E64,E88,E96]' \
    # diffuser.T={diffuser_steps} # DON'T MESS WITH THIS \
    if fix_contigs:
        os.remove(temp_file.name)
    produced_pdbs = [output_dir / f"{output_prefix}_{i}.pdb" for i in range(num_designs)]
    if all(p.exists() for p in produced_pdbs):
        print("All output PDBs already exist, skipping RFDiffusion run.")
        return rfdiffusion_command
    if run:
        subprocess.run(rfdiffusion_command, shell=True, check=True)
        return rfdiffusion_command
    else:
        return rfdiffusion_command


# def run_rfdiffusion_singularity(
#     input_pdb_file: str | Path,
#     output_dir: str | Path,
#     output_prefix: str = "rfdiff_design",
#     singularity_exec_file: str | Path = config.RFDIFFUSION_SINGULARITY,
#     num_designs: int = 8,
#     contig_str: str = "",
#     # diffuser_steps: int = 50,
#     hotspot_residues: str = "",
#     write_trajectory: bool = False,
#     extra_args: str | Path = "",
#     fix_contigs: bool = True,
#     run: bool = True,
# ) -> str:
#     '''
#     # X1-20/2-5/Y76-94/0 B1-107
#     '''
#     output_dir = Path(output_dir).resolve()
#     output_dir.mkdir(parents=True, exist_ok=True)
#     if fix_contigs:
#         temp_file, contig_str = fix_contigs_avoid_missing_residues(contig_str, input_pdb_file)
#         contig_str = fix_contigs_avoid_missing_residues(contig_str, temp_file)
#         input_pdb_file = Path(temp_file.name).resolve()
#         if (output_dir / input_pdb_file.name).exists():
#             print(f"File {output_dir / input_pdb_file.name} already exists, skipping copy of preprocessed pdb.")
#         else:
#             shutil.copy(input_pdb_file, output_dir / "input.pdb")
#     else:
#         input_pdb_file = Path(input_pdb_file).resolve()
#     input_pdb_dir = input_pdb_file.parent
#     input_pdb_filename = input_pdb_file.name
#     schedule_path = output_dir / "schedules"
#     schedule_path.mkdir(parents=True, exist_ok=True)
#     rfdiffusion_command = f"""singularity exec --nv \
#     --bind {output_dir}:/container/output \
#     --bind {schedule_path}:/container/schedules \
#     --bind {input_pdb_dir}:/container/input \
#     {singularity_exec_file} \
#     python /app/RFdiffusion/scripts/run_inference.py \
#     inference.output_prefix=/container/output/{output_prefix} \
#     inference.model_directory_path=/app/RFdiffusion/models \
#     inference.schedule_directory_path=/container/schedules \
#     inference.input_pdb=/container/input/{input_pdb_filename} \
#     inference.num_designs={num_designs} \
#     'contigmap.contigs=["{contig_str}"]' \
#     inference.write_trajectory={str(write_trajectory).lower()} \
#     {extra_args}"""
#     # ppi.hotspot_residues='[E64,E88,E96]' \
#     # diffuser.T={diffuser_steps} # DON'T MESS WITH THIS \
#     produced_pdbs = [output_dir / f"{output_prefix}_{i}.pdb" for i in range(num_designs)]
#     if all(p.exists() for p in produced_pdbs):
#         print("All output PDBs already exist, skipping RFDiffusion run.")
#         return rfdiffusion_command
#     if run:
#         subprocess.run(rfdiffusion_command, shell=True, check=True)
#         contigs = contig_str.replace('/0 ', ' ').split(" ")
#         # contigs = contig_str.split(" ")
#         for pdb in produced_pdbs:
#             with open(pdb, "r") as f:
#                 pdb_str = f.read()
#             fixed_pdb_str = colabdesign_utils.fix_pdb(pdb_str, contigs)
#             with open(pdb, "w") as f:
#                 f.write(fixed_pdb_str)
#         return rfdiffusion_command
#     else:
#         return rfdiffusion_command

