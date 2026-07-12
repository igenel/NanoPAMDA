#!/usr/bin/env python3
"""
Retain only Cas protein hits located within +/-20 kb of a detected CRISPR
array on the same contig (Methods 1.2 genomic proximity criterion).

Usage:
    python 09_filter_20kb_proximity.py <contig_prot_dir> <piler_summary_dir> <output_dir>
"""
import os
import re
import sys


def process_contig_prot(contig_prot_directory, piler_directory, output_directory):
    os.makedirs(output_directory, exist_ok=True)

    contig_prot_files = [f for f in os.listdir(contig_prot_directory) if f.endswith('_contig_prot')]

    for contig_prot_filename in contig_prot_files:
        contig_prot_file_path = os.path.join(contig_prot_directory, contig_prot_filename)
        output_file_path = os.path.join(output_directory, contig_prot_filename.replace('_contig_prot', '_output.txt'))

        piler_filename = re.sub(r'.*_hmm_results_', '', contig_prot_filename).replace('_contig_prot', '_pilercr.out.sum.start_end')
        piler_file_path = os.path.join(piler_directory, piler_filename)

        if not os.path.exists(piler_file_path):
            print(f"No corresponding piler file found for {contig_prot_filename}")
            continue

        piler_data = {}
        with open(piler_file_path, 'r') as piler_file:
            for line in piler_file:
                columns = line.strip().split()
                piler_data[columns[0]] = (int(columns[1]), int(columns[2]))

        with open(contig_prot_file_path, 'r') as contig_prot_file, open(output_file_path, 'w') as output_file:
            for line in contig_prot_file:
                columns = line.strip().split()
                contig_name = columns[0]
                protein_start = int(columns[2])
                protein_end = int(columns[3])

                if contig_name in piler_data:
                    array_start, array_end = piler_data[contig_name]
                    if abs(protein_start - array_start) <= 20000 and abs(protein_end - array_end) <= 20000:
                        output_file.write(line)


if __name__ == "__main__":
    contig_prot_directory = sys.argv[1]
    piler_directory = sys.argv[2]
    output_directory = sys.argv[3]
    process_contig_prot(contig_prot_directory, piler_directory, output_directory)
