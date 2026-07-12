#!/bin/bash
# Extract CRISPR-array-containing contigs and predict ORFs/proteins on them
# using Prodigal v2.6.3 in metagenomic mode (Methods 1.2).
#
# Requires `prodigal` to be available on $PATH (e.g., via a conda environment
# activated before running this script).
#
# Usage: ./03_run_prodigal.sh <piler_summary_dir> <renamed_contigs_dir> <output_dir>
#   piler_summary_dir : output of 02_run_pilercr.sh, e.g. <piler_out>/summary/loc_start_end
#   renamed_contigs_dir : directory of *.fna files produced by 01_rename_headers.py

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <summary_dir> <input_dir> <output_dir>"
    exit 1
fi

summary_dir="$1"
input_dir="$2"
output_dir="$3"

mkdir -p "$output_dir"

for filepath in "$summary_dir"/*.start_end; do
    if [ -s "$filepath" ]; then
        filename=$(basename "$filepath")
        file_name="${filename%_pilercr.out.sum.start_end}"
        fna_file_path="$input_dir/$file_name.fna"

        contig_ids=()
        while IFS= read -r line; do
            contig_id=$(echo "$line" | awk '{print $1}')
            contig_ids+=("$contig_id")
        done < "$filepath"

        # Extract only the CRISPR-array-containing contigs before running Prodigal
        temp_fna_path="$output_dir/${file_name}_contigs.fna"
        awk -v contigs="$(IFS=,; echo "${contig_ids[*]}")" '
            BEGIN { split(contigs, a, ","); for (i in a) ids[a[i]] }
            /^>/ { id=substr($1, 2); capture = id in ids }
            capture { print }
        ' "$fna_file_path" > "$temp_fna_path"

        gbk_file="$output_dir/$file_name.gbk"
        prot_file="$output_dir/$file_name.faa"

        prodigal -i "$temp_fna_path" -o "$gbk_file" -a "$prot_file" -p meta

        rm "$temp_fna_path"
    fi
done
