#!/usr/bin/env python3
"""
Merge CD-HIT cluster assignments back into the main candidate metadata table
(Methods 1.4).

Usage:
    python 04_merge_cluster_assignments.py <candidates.tsv> <cluster_map.tsv> <output.tsv>

Both input files must share a 'unique_id' column to join on.
"""
import sys
import pandas as pd

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python 04_merge_cluster_assignments.py <candidates.tsv> <cluster_map.tsv> <output.tsv>")
        sys.exit(1)

    candidates_file = sys.argv[1]
    cluster_map_file = sys.argv[2]
    output_file = sys.argv[3]

    df1 = pd.read_csv(candidates_file, sep='\t')
    df2 = pd.read_csv(cluster_map_file, sep='\t')

    merged_df = pd.merge(df1, df2, on='unique_id', how='inner')
    merged_df.to_csv(output_file, sep='\t', index=False)

    print(f"Files merged successfully! -> {output_file}")
