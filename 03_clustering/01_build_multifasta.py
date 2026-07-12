#!/usr/bin/env python3
"""
Convert a length-filtered, unique_id-assigned Cas9 candidate TSV into
multiFASTA format for CD-HIT clustering (Methods 1.4).

Each record's header is its 'unique_id' column value (e.g. 'seq_1'), so that
downstream cluster-map extraction (03_extract_cluster_map.sh) can correctly
parse CD-HIT's .clstr output back to the candidate table.

Usage:
    python 01_build_multifasta.py <input.tsv> <output.fa>

Expects the input TSV to have 'unique_id' and 'sequence' columns (output of
00_assign_unique_ids.py).
"""
import csv
import sys

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python 01_build_multifasta.py <input.tsv> <output.fa>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    with open(input_file, 'r') as tsvfile, open(output_file, 'w') as fastafile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for row in reader:
            unique_id = row['unique_id']
            sequence = row['sequence']
            fastafile.write(f">{unique_id}\n{sequence}\n")

    print(f"Multifasta file '{output_file}' generated successfully.")
