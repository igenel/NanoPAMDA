#!/bin/bash
# Join merged HMM hits with their CDS genomic coordinates on the shared
# protein ID field (Methods 1.2), producing one contig/protein/location
# table per dataset.
#
# Usage: ./08_join_contig_protein.sh <merged_hits_dir> <protein_locations_dir> <output_dir>

set -euo pipefail

merged_dir=$1
locations_dir=$2
output_dir=$3

mkdir -p "$output_dir"

for merged_file in "${merged_dir}"/*_merged.out; do
    base_name=$(basename "${merged_file}" _merged.out)

    merged_file_path="${merged_dir}/${base_name}_merged.out"
    locations_file_path="${locations_dir}/${base_name}_merged_locations_output.txt"

    if [[ ! -f "$merged_file_path" ]]; then
        echo "File not found: $merged_file_path"
        continue
    fi
    if [[ ! -f "$locations_file_path" ]]; then
        echo "File not found: $locations_file_path"
        continue
    fi

    output_file="${output_dir}/${base_name}_contig_prot"

    sort -k2,2 "$merged_file_path" > "${merged_file_path}.sorted"
    sort -k1,1 "$locations_file_path" > "${locations_file_path}.sorted"

    join -1 2 -2 1 -o 1.1,1.2,2.2,2.3 "${merged_file_path}.sorted" "${locations_file_path}.sorted" > "$output_file"

    rm "${merged_file_path}.sorted" "${locations_file_path}.sorted"
done
