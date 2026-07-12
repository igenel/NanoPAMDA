#!/usr/bin/env python3
"""
Assign a simple sequential unique_id (seq_1, seq_2, ...) to every candidate
in the collected candidate table, once all data (MAG-derived + dbGaP-derived,
post length-filtering) has been merged (Methods 1.4).

This unique_id is what downstream clustering (CD-HIT), the cluster map
extraction, and the merge-back step all key on -- it is independent of the
'natural_<i>_<cas_class>' style label used elsewhere for readability.

Usage:
    python 00_assign_unique_ids.py <input.tsv> <output.tsv>

Adds a 'unique_id' column as the first column of the output TSV.
"""
import sys
import pandas as pd

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python 00_assign_unique_ids.py <input.tsv> <output.tsv>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    df = pd.read_csv(input_file, sep='\t')
    df.insert(0, 'unique_id', [f"seq_{i}" for i in range(1, len(df) + 1)])
    df.to_csv(output_file, sep='\t', index=False)

    print(f"Assigned {len(df)} unique_id values (seq_1 .. seq_{len(df)}) -> {output_file}")
