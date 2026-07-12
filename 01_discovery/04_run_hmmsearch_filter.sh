#!/bin/bash
# Screen predicted proteins against the curated Cas profile HMM library using
# hmmsearch (HMMER v3.3.2), then apply the alignment coverage/bitscore/e-value
# retention thresholds described in Methods 1.2 (coverage >70%, bitscore >=40,
# e-value <1e-5).
#
# Requires `hmmsearch` on $PATH.
#
# Usage: ./04_run_hmmsearch_filter.sh <hmm_profiles_dir> <prodigal_results_dir> <output_dir>
#   hmm_profiles_dir: directory containing subdirectories of *.hmm profile files
#   prodigal_results_dir: directory of *.faa files from 03_run_prodigal.sh
#
# Note: profiles are processed in small batches (default 4 subdirectories at a
# time) purely to make progress easier to monitor on large runs; adjust
# BATCH_SIZE below or set the environment variable of the same name.

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <profiles_dir> <prodigal_results_dir> <output_dir>"
    exit 1
fi

profiles_dir=$1
prodigal_results_dir=$2
output_dir=$3
BATCH_SIZE="${BATCH_SIZE:-4}"
HMMSEARCH_CPU="${HMMSEARCH_CPU:-4}"

mkdir -p "$output_dir"

subdirs=($(ls -d "$profiles_dir"/* | xargs -n 1 basename))
num_subdirs=${#subdirs[@]}

for (( start_index = 0; start_index < num_subdirs; start_index += BATCH_SIZE )); do
    end_index=$(( start_index + BATCH_SIZE - 1 ))
    if [ $end_index -ge $num_subdirs ]; then
        end_index=$(( num_subdirs - 1 ))
    fi

    echo "Running batch from index $start_index to $end_index"

    for (( i = start_index; i <= end_index; i++ )); do
        subdir_name=${subdirs[$i]}
        echo "Processing subdirectory: $subdir_name"
        subdir="$profiles_dir/$subdir_name"

        if [ -d "$subdir" ]; then
            output_subdir="$output_dir/${subdir_name}_hmm_results"
            mkdir -p "$output_subdir"

            hmm_files=("$subdir"/*.hmm)
            if [ ${#hmm_files[@]} -eq 0 ]; then
                echo "No .hmm files found in subdirectory $subdir."
                continue
            fi

            for hmm_file in "${hmm_files[@]}"; do
                hmm_name=$(basename "$hmm_file" .hmm)
                for faa_file in "$prodigal_results_dir"/*.faa; do
                    faa_name=$(basename "$faa_file" .faa)
                    output_domtbl="$output_subdir/${hmm_name}_${faa_name}.domtblout"
                    echo "Running hmmsearch for $hmm_file against $faa_file..."
                    hmmsearch --cpu "$HMMSEARCH_CPU" --domtblout "$output_domtbl" "$hmm_file" "$faa_file" > /dev/null

                    if [ -f "$output_domtbl" ]; then
                        echo "Filtering $output_domtbl..."
                        # Retention criteria (Methods 1.2): alignment coverage >70%,
                        # bitscore >=40, e-value <1e-5
                        awk 'BEGIN { FS=" "; OFS="\t" }
                            $1 ~ /^#/ { next }
                            {
                                evalue = $7;
                                bitscore = $8;
                                qlen = $6;
                                aln_start = $18;
                                aln_end = $19;
                                aln_coverage = ($19 - $18 + 1) / qlen;
                                if (bitscore >= 40 && evalue < 1e-5 && aln_coverage > 0.7) {
                                    print
                                }
                            }' "$output_domtbl" > "$output_subdir/${hmm_name}_${faa_name}_filtered.domtblout"
                        echo "Processed and filtered $output_domtbl"
                    else
                        echo "Output domtblout file $output_domtbl not found."
                    fi
                done
            done
        else
            echo "Subdirectory $subdir not found."
        fi
        echo "Finished processing subdirectory: $subdir_name"
    done
done

echo "All batches processed."
