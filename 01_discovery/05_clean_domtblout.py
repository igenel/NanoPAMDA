#!/usr/bin/env python3
"""
Parse filtered HMMER domtblout files into a compact "contig protein-index"
format for downstream merging (Methods 1.2).

For each *_filtered.domtblout file (output of 04_run_hmmsearch_filter.sh),
extracts the contig index and protein ID pair from each retained hit line.

Usage:
    python 05_clean_domtblout.py <hmm_results_dir> <prodigal_faa_dir> <output_dir>
"""
import os
import re
import sys
from pathlib import Path

input_dir = sys.argv[1]
faa_dir = sys.argv[2]
output_root_dir = sys.argv[3]
Path(output_root_dir).mkdir(parents=True, exist_ok=True)

pattern = r'contig-100_(\d+)_(\d+).*ID=(\d+)_.*;'

file_contents = {}

for root, dirs, files in os.walk(input_dir):
    for filename in files:
        if filename.endswith("_filtered.domtblout"):
            input_path = os.path.join(root, filename)

            relative_dir = os.path.relpath(root, input_dir)
            output_dir = os.path.join(output_root_dir, relative_dir)
            Path(output_dir).mkdir(parents=True, exist_ok=True)

            output_path = os.path.join(output_dir, filename)

            with open(input_path, 'r') as input_file, open(output_path, 'w') as output_file:
                for line in input_file:
                    match = re.search(pattern, line)
                    if match:
                        k_value_1 = match.group(1)
                        id_value_1 = match.group(3)
                        id_value_2 = match.group(2)
                        output_line = f"contig-100_{k_value_1} {id_value_1}_{id_value_2}\n"
                        output_file.write(output_line)

                        if filename in file_contents:
                            file_contents[filename].add(output_line.strip())
                        else:
                            file_contents[filename] = {output_line.strip()}
