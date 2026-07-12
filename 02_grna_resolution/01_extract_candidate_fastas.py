#!/usr/bin/env python3
"""
Extract per-candidate FASTA sequences and build a processing manifest for
CRISPRtracrRNA (Methods 1.3).

For each 20kb-proximity-filtered Cas9 hit (Section 1.2, step 9 output), looks
up the corresponding contig sequence in the renamed source .fna file and
writes it to its own FASTA file, plus records a manifest row (source HMM
file, output FASTA path, contig ID, protein ID, hit coordinates) to a TSV.

Usage:
    python 01_extract_candidate_fastas.py

Adjust base_parent_dir / output paths below (or refactor to accept CLI args)
for your own environment. The original internal version iterated over
subdir_1 .. subdir_30 (one per chunk of a specific dataset); generalize the
loop range/paths as needed for other datasets.
"""
import os
import glob
import pandas as pd
from Bio import SeqIO

# ---------------------------------------------------------------------------
# CONFIG - adjust for your environment / dataset
# ---------------------------------------------------------------------------
base_parent_dir = "./mag_pipeline_150k"  # root containing subdir_1 ... subdir_N
output_dir = os.path.join(base_parent_dir, "tracr_input")
num_subdirs = 30
# ---------------------------------------------------------------------------

os.makedirs(output_dir, exist_ok=True)

tsv_data = []

for i in range(1, num_subdirs + 1):
    hmm_dir = os.path.join(base_parent_dir, f"subdir_{i}", "hmm_results_output", "merged", "prot_loc", "contig_prot", "new_20kb")
    fna_dir = os.path.join(base_parent_dir, f"subdir_{i}")

    for hmm_file in glob.glob(os.path.join(hmm_dir, "cas9_hmm_results_*_output.txt")):
        base_name = os.path.basename(hmm_file).replace("cas9_hmm_results_", "").replace("_output.txt", "")
        fna_file = os.path.join(fna_dir, f"{base_name}.fna")

        if os.path.isfile(fna_file):
            print(f"Processing {hmm_file} and {fna_file}")
            sequences = SeqIO.to_dict(SeqIO.parse(fna_file, "fasta"))

            with open(hmm_file, 'r') as hmm_f:
                for line in hmm_f:
                    contig_id, prot_id, start_num, end_num = line.strip().split()

                    if contig_id in sequences:
                        sequence = str(sequences[contig_id].seq)
                        output_fasta = os.path.join(output_dir, f"{base_name}_{contig_id}.fasta")

                        with open(output_fasta, 'w') as fasta_f:
                            fasta_f.write(f">{base_name}_{contig_id}\n{sequence}\n")

                        tsv_data.append([hmm_file, output_fasta, contig_id, prot_id, start_num, end_num])
                    else:
                        print(f"WARNING: contig_id {contig_id} not found in {fna_file}")
        else:
            print(f"WARNING: .fna file not found for {base_name}")

tsv_file = os.path.join(output_dir, "results.tsv")
df = pd.DataFrame(tsv_data, columns=["HMM_File", "FASTA_File", "Contig_ID", "Protid", "Start_Number", "End_Number"])
df.to_csv(tsv_file, sep='\t', index=False)
print(f"Manifest written: {tsv_file}")
