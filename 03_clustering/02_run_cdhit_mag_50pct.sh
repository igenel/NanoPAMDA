#!/bin/bash
# CD-HIT clustering of MAG-derived candidates at 50% sequence identity
# (Methods 1.4).
#
# Requires `cd-hit` on $PATH.
#
# Usage: ./02_run_cdhit_mag_50pct.sh <input.fa> <output_prefix>

set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.fa> <output_prefix>"
    exit 1
fi

input_fa="$1"
output_prefix="$2"

cd-hit -i "$input_fa" -o "${output_prefix}.fasta" -c 0.50 -n 5

echo "CD-HIT complete: ${output_prefix}.fasta (representatives), ${output_prefix}.fasta.clstr (cluster assignments)"
