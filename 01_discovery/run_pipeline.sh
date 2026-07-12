#!/bin/bash
# Orchestrator for the full CRISPR array / Cas9 candidate discovery pipeline
# (Methods 1.1-1.2). Runs steps 01-10 in order against a single input dataset
# directory of raw contigs.
#
# All hardcoded cluster-specific paths (SLURM account, conda install path,
# personal directories) have been removed/genericized relative to the
# original internal version of this script. Adjust the CONFIG section below
# for your own environment, or export the corresponding environment
# variables before running.
#
# Usage: ./run_pipeline.sh <input_dir>
#   <input_dir> should contain the raw contig FASTA files for one dataset
#   (e.g., one MAG catalog or one dbGaP-derived MetaHipMer2 assembly output).

set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_dir>"
    exit 1
fi

# ---------------------------------------------------------------------------
# CONFIG - adjust for your environment
# ---------------------------------------------------------------------------
input_dir="$1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Path to the directory of Cas profile HMMs (see Methods 1.2 / Makarova et al., 2015)
PROFILES_DIR="${PROFILES_DIR:-./profiles_combined}"

# Optional: activate your own environment(s) providing pilercr, prodigal,
# hmmsearch before running this script, e.g.:
#   conda activate cas9-discovery
# ---------------------------------------------------------------------------

output_dir_piler="$input_dir/piler_output"
summary_dir="$output_dir_piler/summary/loc_start_end"
output_dir_prodigal="$input_dir/prodigal_output"
output_dir_hmm="$input_dir/hmm_results"
output_root_dir="$input_dir/hmm_results_output"
merged_out_dir="$output_root_dir/merged"
prot_loc_dir="$output_root_dir/merged/prot_loc"
contig_prot_dir="$output_root_dir/merged/prot_loc/contig_prot"
new_20kb_dir="$output_root_dir/merged/prot_loc/contig_prot/new_20kb"
tsv_dir="$input_dir/multifasta"

echo "Step 1/10: Standardizing contig headers..."
python "$SCRIPT_DIR/01_rename_headers.py" "$input_dir"

echo "Step 2/10: Running PILER-CR..."
bash "$SCRIPT_DIR/02_run_pilercr.sh" "$input_dir" "$output_dir_piler"

echo "Step 3/10: Running Prodigal..."
bash "$SCRIPT_DIR/03_run_prodigal.sh" "$summary_dir" "$input_dir" "$output_dir_prodigal"

echo "Step 4/10: Running hmmsearch + filtering..."
bash "$SCRIPT_DIR/04_run_hmmsearch_filter.sh" "$PROFILES_DIR" "$output_dir_prodigal" "$output_dir_hmm"

echo "Step 5/10: Cleaning domtblout hits..."
python "$SCRIPT_DIR/05_clean_domtblout.py" "$output_dir_hmm" "$output_dir_prodigal" "$output_root_dir"

echo "Step 6/10: Merging domtblout hits..."
python "$SCRIPT_DIR/06_merge_domtblout.py" "$output_root_dir"

echo "Step 7/10: Extracting protein CDS locations..."
python "$SCRIPT_DIR/07_get_protein_locations.py" "$merged_out_dir" "$output_dir_prodigal" "$prot_loc_dir"

echo "Step 8/10: Joining contigs with protein locations..."
bash "$SCRIPT_DIR/08_join_contig_protein.sh" "$merged_out_dir" "$prot_loc_dir" "$contig_prot_dir"

echo "Step 9/10: Filtering to +/-20kb of CRISPR arrays..."
python "$SCRIPT_DIR/09_filter_20kb_proximity.py" "$contig_prot_dir" "$summary_dir" "$new_20kb_dir"

echo "Step 10/10: Building final Cas9 candidate MultiFASTA/TSV database..."
python "$SCRIPT_DIR/10_build_cas_tsv.py" "$new_20kb_dir" "$output_dir_prodigal" "$tsv_dir"

echo "Pipeline complete. Final candidate database: $tsv_dir"
