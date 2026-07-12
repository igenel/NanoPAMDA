#!/bin/bash
# Extract a simple two-column (unique_id, cluster_id) map from a CD-HIT
# .clstr file (Methods 1.4).
#
# Usage: ./03_extract_cluster_map.sh <clusters.fasta.clstr> <output_map.tsv> <cluster_id_column_name>
#   cluster_id_column_name: name for the second column, e.g. "50_id" or "90_id"
#   (matches the identity threshold used to generate the .clstr file, for
#   downstream tracking across the MAG (50%) and dbGaP (90%) arms)

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <clusters.fasta.clstr> <output_map.tsv> <cluster_id_column_name>"
    exit 1
fi

clstr_file="$1"
output_file="$2"
id_column_name="$3"

{
  echo -e "unique_id\t${id_column_name}"
} > "$output_file"

awk -v col="$id_column_name" '
    /^>/ { cluster_id++; next }
    /^[0-9]/ {
        match($0, />seq_[^\.]+/, seq)
        gsub(">", "", seq[0])
        print seq[0] "\t" cluster_id
    }
' "$clstr_file" >> "$output_file"

echo "$output_file created successfully."
