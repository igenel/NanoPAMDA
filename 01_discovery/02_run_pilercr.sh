#!/bin/bash
# Screen assemblies for CRISPR arrays using PILER-CR v1.06 (Methods 1.2).
#
# Requires the `pilercr` binary to be available. By default this script looks
# for it on $PATH; override by setting the PILERCR_BIN environment variable,
# e.g.: export PILERCR_BIN=/path/to/pilercr1.06/pilercr
#
# Usage: ./02_run_pilercr.sh <input_dir_with_fna_files> <output_dir>

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    exit 1
fi

input_dir="$1"
output_dir="$2"
PILERCR_BIN="${PILERCR_BIN:-pilercr}"

mkdir -p "$output_dir"

# Step 1: run PILER-CR on every renamed contig file
for contigs_file in "$input_dir"/*.fna; do
    if [ -f "$contigs_file" ]; then
        filename=$(basename "$contigs_file" .fna)
        m_output_file="$output_dir/${filename}_pilercr.out"
        "$PILERCR_BIN" -in "$contigs_file" -out "$m_output_file"
    else
        echo "Error: $contigs_file not found in $input_dir"
    fi
done

# Step 2: extract the "SUMMARY BY SIMILARITY" block from each PILER-CR output
summary_dir="$output_dir/summary"
mkdir -p "$summary_dir"

for input_file in "$output_dir"/*_pilercr.out; do
    output_file="${summary_dir}/${input_file##*/}.sum"
    awk '/SUMMARY BY SIMILARITY/{flag=1; next} /SUMMARY BY POSITION/{flag=0} flag' "$input_file" > "$output_file"
done

# Step 3: extract contig ID + array start/end coordinates for each hit
start_end_dir="$summary_dir/loc_start_end"
mkdir -p "$start_end_dir"

for n_input_file in "$summary_dir"/*.sum; do
    n_output_file="${start_end_dir}/${n_input_file##*/}.start_end"
    awk '
        /^[[:space:]]*[0-9]+/ {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^[0-9]+$/ && $(i+1) ~ /^[0-9]+$/) {
                    start = $i
                    end = $(i+1) + start
                    for (j=1; j<=NF; j++) {
                        if ($j ~ /^contig-100_[0-9]+$/) {
                            contig_part = $j
                            break
                        }
                    }
                    print contig_part, start, end
                    break
                }
            }
        }
    ' "$n_input_file" > "$n_output_file"
done
