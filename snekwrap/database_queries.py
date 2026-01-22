import requests
import os
from pathlib import Path
import pandas as pd


def download_pdb(pdb_id, output_dir: str|Path ='.'):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # Ensure PDB ID is lowercase
    pdb_id = pdb_id.lower()
    
    # Construct the URL
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    
    # Send a GET request to the URL
    response = requests.get(url)
    
    # Check if the request was successful
    if response.status_code == 200:
        # Construct the output file path
        output_file = output_dir / f"{pdb_id}.pdb"
        
        # Write the content to a file
        with open(output_file, 'wb') as f:
            f.write(response.content)
        
        print(f"Successfully downloaded {pdb_id}.pdb")
        return output_file
    else:
        raise ValueError(f"Failed to download {pdb_id}.pdb. Status code: {response.status_code}")


def download_cif(pdb_id, output_dir: str|Path ='.'):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    # Ensure PDB ID is lowercase
    pdb_id = pdb_id.lower()
    
    # Construct the URL
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    
    # Send a GET request to the URL
    response = requests.get(url)
    
    # Check if the request was successful
    if response.status_code == 200:
        # Construct the output file path
        output_file = output_dir / f"{pdb_id}.cif"
        
        # Write the content to a file
        with open(output_file, 'wb') as f:
            f.write(response.content)
        
        print(f"Successfully downloaded {pdb_id}.cif")
        return output_file
    else:
        raise ValueError(f"Failed to download {pdb_id}.cif. Status code: {response.status_code}")


def download_uniprot_sequence(uniprot_id, output_dir: str|Path|None = None):
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / f"{uniprot_id}.fasta"
        if output_file.exists():
            print(f"File {output_file} already exists. Skipping download.")
            return output_file
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        if output_dir is None:
            # Split the response text into lines
            lines = response.text.split('\n')
            # Join all lines except the first one (header) and remove whitespace
            sequence = ''.join(lines[1:]).replace(' ', '')
            return sequence
        with open(output_file, 'w') as f:
            f.write(response.text)
        return output_file
    else:
        print(f"Failed to download sequence for {uniprot_id}. Status code: {response.status_code}")
        return None


def download_uniprot_gff(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.gff"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        raise Exception(f"Failed to download GFF for {uniprot_id}. Status code: {response.status_code}")


def parse_gff_to_dataframe(gff_content):
    # Skip comment lines and split the content into lines
    lines = [line for line in gff_content.split('\n') if line and not line.startswith('#')]
    
    # Define column names for GFF format
    columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    
    # Parse lines into a list of dictionaries
    data = []
    for line in lines:
        fields = line.strip().split('\t')
        assert len(fields) == 9  # Ensure we have all 9 fields
        row = dict(zip(columns, fields))
        # Parse attributes
        attrs = dict(item.split('=') for item in row['attributes'].split(';') if '=' in item)
        row.update(attrs)
        data.append(row)

    # Create DataFrame
    df = pd.DataFrame(data)
    return df


def get_uniprot_annotations_from_gff(uniprot_id):
    gff_content = download_uniprot_gff(uniprot_id)
    df = parse_gff_to_dataframe(gff_content)
    return df