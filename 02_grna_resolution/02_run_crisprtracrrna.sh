#!/bin/bash
# Organize extracted candidate FASTA files into per-file subdirectories and
# run CRISPRtracrRNA (BackofenLab; model_type II) on each (Methods 1.3).
#
# CRISPRtracrRNA (https://github.com/BackofenLab/CRISPRtracrRNA) must be
# installed separately, including its CRISPRidentify and CRISPRcasIdentifier
# dependencies -- see that repository's README for setup instructions and its
# own conda environment ("crispr_tracr_rna_env").
#
# Usage: ./02_run_crisprtracrrna.sh <input_dir_with_fasta_files> <output_dir> <path_to_CRISPRtracrRNA.py>

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_dir> <output_dir> <crisprtracrrna_script_path>"
    exit 1
fi

input_dir="$1"
output_dir="$2"
CRISPRTRACRRNA_PY="$3"

mkdir -p "$output_dir"

# Step 1: distribute each .fasta file into its own subdir_N (CRISPRtracrRNA
# expects one input folder per run)
counter=1
for fasta_file in "$input_dir"/*.fasta; do
    if [[ ! -e "$fasta_file" ]]; then
        echo "No .fasta files found in $input_dir"
        exit 1
    fi
    subdir_name="subdir_$counter"
    subdir_path="$input_dir/$subdir_name"
    mkdir -p "$subdir_path"
    mv "$fasta_file" "$subdir_path/"
    counter=$((counter + 1))
done

echo "Distributed FASTA files into per-candidate subdirectories under $input_dir"

# Step 2: run CRISPRtracrRNA on each subdirectory
for subdir in "$input_dir"/subdir_*; do
    if [[ -d "$subdir" ]]; then
        subdir_name=$(basename "$subdir")
        python "$CRISPRTRACRRNA_PY" --input_folder "$subdir" --output_folder "$output_dir" --model_type II
    else
        echo "Directory $subdir does not exist."
    fi
done

echo "CRISPRtracrRNA processing complete. Output: $output_dir"
