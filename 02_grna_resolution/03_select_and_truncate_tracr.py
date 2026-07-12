#!/usr/bin/env python3
"""
Select one representative tracrRNA candidate per locus from CRISPRtracrRNA
output, apply the terminator-truncation heuristic, and pair each resolved
tracrRNA with its corresponding Cas9 protein sequence (Methods 1.3).

Selection logic per locus:
  1. Among all CRISPRtracrRNA-nominated array elements, retain only the
     element(s) with the greatest crispr_array_start (closest to the
     candidate tracrRNA locus).
  2. Among ties, select the element with the lowest interaction_energy.

Terminator truncation:
  - Search for the first 'TTTT' occurring after nucleotide position 40 of
    the tracrRNA sequence; truncate there if found. No upper bound is
    applied. Sequences without a hit are kept at full length.

Usage:
    python 03_select_and_truncate_tracr.py
(adjust the hardcoded paths in the __main__ block below for your environment,
or refactor to accept CLI arguments)
"""
import os
import csv


def process_csv(csv_path, new_name, protid, cas9_sequences):
    if not os.path.exists(csv_path):
        print(f"File not found: {csv_path}, skipping.")
        return None

    with open(csv_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        rows = list(reader)

        if not rows:
            print(f"No data in file: {csv_path}, skipping.")
            return None

        # Step 1: retain element(s) with greatest crispr_array_start
        max_start = max(int(row["crispr_array_start"]) for row in rows)
        rows_with_max_start = [row for row in rows if int(row["crispr_array_start"]) == max_start]

        # Step 2: among ties, lowest interaction_energy
        min_energy_row = min(
            rows_with_max_start,
            key=lambda row: float(row["interaction_energy"]) if row["interaction_energy"] else float('inf')
        )

        accession_number = min_energy_row["accession_number"]
        crispr_array_repeat_consensus = min_energy_row["crispr_array_repeat_consensus"]
        crispr_array_orientation = min_energy_row["crispr_array_orientation"]
        anti_repeat_sequence = min_energy_row["anti_repeat_sequence"]
        tracr_rna_sequence = min_energy_row["tracr_rna_sequence"]

        final_tracr, truncated_flag = generate_final_tracr(tracr_rna_sequence)

        cas9_sequence = cas9_sequences.get((new_name, protid), None)

        return {
            "accession_number": accession_number,
            "crispr_array_repeat_consensus": crispr_array_repeat_consensus,
            "crispr_array_orientation": crispr_array_orientation,
            "anti_repeat_sequence": anti_repeat_sequence,
            "tracr_rna_sequence": tracr_rna_sequence,
            "final_tracr": final_tracr,
            "truncated_flag": truncated_flag,
            "cas9_sequence": cas9_sequence,
            "protid": protid
        }


def generate_final_tracr(tracr_rna_sequence):
    """
    Truncate at the first 'TTTT' occurring after nucleotide position 40.
    No upper bound is applied (Arnold was not used for this step).
    """
    index = tracr_rna_sequence.find('TTTT', 40)
    if index != -1:
        final_tracr = tracr_rna_sequence[:index]
        truncated_flag = "Yes"
    else:
        final_tracr = tracr_rna_sequence
        truncated_flag = "No"
    return final_tracr, truncated_flag


def load_cas9_sequences(merged_output_file):
    """
    Load Cas9 protein sequences from the Section 1.2 final candidate TSV
    (output of 01_discovery/10_build_cas_tsv.py), keyed by
    (source_file_base_name, protein_id).
    """
    cas9_sequences = {}
    with open(merged_output_file, newline='') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for row in reader:
            if row["cas_class"] == "cas9":
                new_name = row["source_file_base_name"].split("_contig")[0]
                protid = row["protein_id"]
                sequence = row["sequence"]
                cas9_sequences[(new_name, protid)] = sequence
    return cas9_sequences


if __name__ == "__main__":
    # ---- CONFIG - adjust for your environment ----
    results_file = "./tracr_input/results.tsv"                 # from 01_extract_candidate_fastas.py
    output_dir = "./tracr_output/"                              # CRISPRtracrRNA per-candidate CSV outputs
    output_file = "./tracr_processed_output.tsv"
    merged_output_file = "./merged_output_files/all_sequences.tsv"  # from 01_discovery/10_build_cas_tsv.py
    # -----------------------------------------------

    cas9_sequences = load_cas9_sequences(merged_output_file)

    with open(results_file, newline='') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')

        with open(output_file, mode='w', newline='') as outfile:
            fieldnames = [
                "accession_number", "crispr_array_repeat_consensus",
                "crispr_array_orientation", "anti_repeat_sequence",
                "tracr_rna_sequence", "final_tracr", "truncated_flag",
                "cas9_sequence", "protid"
            ]
            writer = csv.DictWriter(outfile, fieldnames=fieldnames)
            writer.writeheader()

            for row in reader:
                fasta_file_path = row["FASTA_File"]
                file_name = os.path.basename(fasta_file_path).replace(".fasta", "")
                new_name = file_name.split("_contig")[0]
                protid = row["Protid"]
                csv_path = os.path.join(output_dir, f"{file_name}.csv")

                result = process_csv(csv_path, new_name, protid, cas9_sequences)
                if result:
                    writer.writerow(result)
                else:
                    print(f"Skipped {csv_path}")
