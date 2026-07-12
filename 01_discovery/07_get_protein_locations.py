#!/usr/bin/env python3
"""
Extract CDS genomic coordinates for each retained Cas protein hit from the
corresponding Prodigal-generated GenBank (.gbk) file (Methods 1.2).

Usage:
    python 07_get_protein_locations.py <merged_hits_dir> <prodigal_gbk_dir> <output_dir>
"""
import os
import re
import sys


def read_protein_ids(file_path):
    protein_ids = []
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            protein_id = parts[1]
            protein_ids.append(protein_id)
    return protein_ids


def extract_cds_locations(genbank_file_path, protein_ids):
    cds_locations = {}
    with open(genbank_file_path, 'r') as f:
        content = f.read()
        for protein_id in protein_ids:
            regex_id = fr'(\b{re.escape(protein_id)}\b)'
            match_id = re.search(regex_id, content)
            if match_id:
                start_idx = match_id.start()
                pattern = r'(complement\()?<?(\d+)\.\.>?(\d+)\)?'
                search_start = max(start_idx - 100, 0)
                substring_to_search = content[search_start:start_idx]
                match = re.search(pattern, substring_to_search)
                if match:
                    x, y = match.group(2), match.group(3)
                    cds_locations.setdefault(protein_id, []).append((x, y))
    return cds_locations


def write_output(output_file, locations):
    output_dir = os.path.dirname(output_file)
    os.makedirs(output_dir, exist_ok=True)
    with open(output_file, 'w') as f:
        for protein_id, location_list in locations.items():
            for start, end in location_list:
                f.write(f"{protein_id} {start} {end}\n")


def find_matching_gbk_file(protein_ids_file, genbank_directory):
    base_name = os.path.basename(protein_ids_file)
    name_parts = base_name.split("_")
    file_prefix = name_parts[-2]
    for file in os.listdir(genbank_directory):
        if file.startswith(file_prefix) and file.endswith(".gbk"):
            return os.path.join(genbank_directory, file)
    return None


if __name__ == "__main__":
    protein_ids_directory = sys.argv[1]
    genbank_directory = sys.argv[2]
    output_directory = sys.argv[3]
    os.makedirs(output_directory, exist_ok=True)

    for file in os.listdir(protein_ids_directory):
        if file.endswith("_merged.out"):
            protein_ids_file = os.path.join(protein_ids_directory, file)
            genbank_file = find_matching_gbk_file(protein_ids_file, genbank_directory)
            if not genbank_file:
                print(f"Matching GenBank file not found for '{protein_ids_file}'. Skipping...")
                continue

            output_file = os.path.join(output_directory, f"{os.path.splitext(file)[0]}_locations_output.txt")

            protein_ids = read_protein_ids(protein_ids_file)
            cds_locations = extract_cds_locations(genbank_file, protein_ids)

            if cds_locations:
                write_output(output_file, cds_locations)
                print(f"Output file '{output_file}' successfully created.")
            else:
                print(f"No CDS locations found for '{protein_ids_file}'")
