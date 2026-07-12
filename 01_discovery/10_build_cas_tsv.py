#!/usr/bin/env python3
"""
Build the final centralized MultiFASTA and TSV database of intact
(partial=00) Cas9 candidates (Methods 1.2, final step).

Takes the 20kb-proximity-filtered hit list (output of 09_filter_20kb_proximity.py,
files named *_output.txt) together with the Prodigal .faa files, and:
  1. Parses each hit file to recover the Cas class/subtype (from the parent
     HMM profile identifier embedded in the filename) and the associated
     protein IDs.
  2. Extracts the corresponding full-length protein sequences from the
     Prodigal .faa files, retaining only intact, non-partial ORFs
     (Prodigal partial=00 flag).
  3. Writes all retained sequences to a single MultiFASTA file and a
     matching TSV (cas_class, protein_id, sequence, source_file_base_name).

Usage:
    python 10_build_cas_tsv.py <20kb_filtered_hits_dir> <prodigal_faa_dir> <output_dir>
"""
import os
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import csv


def extract_sequences(faa_file_path, protein_ids, file_base):
    sequences = []
    for record in SeqIO.parse(faa_file_path, "fasta"):
        for protein_id in protein_ids:
            if f"ID={protein_id};partial=00" in record.description:
                # Remove trailing asterisks
                record.seq = record.seq.rstrip("*")
                sequences.append((file_base, protein_id, str(record.seq)))
    return sequences


def process_hmm_results(hmm_results_dir, prodigal_results_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    main_dir_name = os.path.basename(os.path.dirname(output_dir))

    # Step 1: Identify file names and extract protein names, grouped by Cas class
    protein_data = {}

    for file_name in os.listdir(hmm_results_dir):
        if file_name.endswith("_output.txt"):
            # Cas class/subtype is derived from the parent HMM profile identifier
            # embedded in the filename (Methods 1.2)
            cas_class = file_name.split("_hmm_results_")[0].lower().replace('_', '')
            file_name_base = file_name.split("_hmm_results_")[1].split("_output.txt")[0]

            with open(os.path.join(hmm_results_dir, file_name)) as file:
                protein_names = [line.split()[1] for line in file]
                if cas_class not in protein_data:
                    protein_data[cas_class] = {}
                protein_data[cas_class][file_name_base] = protein_names

    print("Extracted protein data:")
    for cas_class, files in protein_data.items():
        for file_base, protein_ids in files.items():
            print(f"cas_class: {cas_class}, file_base: {file_base}, protein_ids: {protein_ids}")

    # Step 2: Collect intact-ORF sequences for each hit
    all_sequences = []

    for cas_class, files in protein_data.items():
        for file_base, protein_ids in files.items():
            faa_file_path = os.path.join(prodigal_results_dir, f"{file_base}.faa")

            print(f"Processing file: {faa_file_path}")
            if os.path.exists(faa_file_path):
                sequences = extract_sequences(faa_file_path, protein_ids, file_base)
                if sequences:
                    print(f"Found sequences for {file_base}: {[seq[1] for seq in sequences]}")
                else:
                    print(f"No sequences found for {file_base}")
                all_sequences.extend([(cas_class, protein_id, seq, file_base) for file_base, protein_id, seq in sequences])
            else:
                print(f"File not found: {faa_file_path}")

    # Step 3: Write MultiFASTA
    output_file_name = f"{main_dir_name}_all_sequences.fasta"
    output_file_path = os.path.join(output_dir, output_file_name)

    seq_records = [SeqRecord(Seq(seq), id=protein_id, description="") for _, protein_id, seq, _ in all_sequences]

    with open(output_file_path, "w") as output_file:
        SeqIO.write(seq_records, output_file, "fasta")

    print(f"Merged MultiFASTA file '{output_file_path}' created successfully.")

    # Step 4: Write matching TSV
    tsv_file_name = f"{main_dir_name}_sequences.tsv"
    tsv_file_path = os.path.join(output_dir, tsv_file_name)

    with open(tsv_file_path, "w", newline="") as tsv_file:
        tsv_writer = csv.writer(tsv_file, delimiter='\t')
        tsv_writer.writerow(['cas_class', 'protein_id', 'sequence', 'source_file_base_name'])
        for cas_class, protein_id, sequence, source_file_base_name in all_sequences:
            tsv_writer.writerow([cas_class, protein_id, sequence, source_file_base_name])

    print(f"TSV file '{tsv_file_path}' created successfully.")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python 10_build_cas_tsv.py <hmm_results_dir> <prodigal_results_dir> <output_dir>")
        sys.exit(1)

    hmm_results_dir = sys.argv[1]
    prodigal_results_dir = sys.argv[2]
    output_dir = sys.argv[3]

    process_hmm_results(hmm_results_dir, prodigal_results_dir, output_dir)
