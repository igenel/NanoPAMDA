#!/usr/bin/env python3
"""
Orchestrator for the full ONT-adapted HT-PAMDA analysis pipeline
(Methods 1.12-1.15): PAM extraction -> baseline correction ->
kinetic modeling -> heatmap generation.

Each step mirrors the corresponding function in the official HT-PAMDA
pipeline (Walton et al., 2021) but is reimplemented here to work on
demultiplexed single-end ONT FASTQs (output of 00_demux_htpamda.py)
rather than paired-end Illumina reads, using a separately-processed
Illumina untreated library as the t=0 baseline control.

Usage:
    python run_pipeline.py

All hardcoded cluster-specific paths from the original internal version
have been replaced with placeholders in the CONFIG section below -- edit
these for each new sequencing run.
"""
import os
import shutil
import pandas as pd

from importlib import import_module

extract_pam_counts = import_module("01_extract_pam_counts")
baseline_correction = import_module("02_baseline_correction")
kinetic_modeling = import_module("03_kinetic_modeling")
generate_heatmaps = import_module("04_generate_heatmaps")

# ══════════════════════════════════════════════════════════════════════════
# CONFIGURATION - edit these for each run
# ══════════════════════════════════════════════════════════════════════════

RUN_NAME = 'ONT_run_name_here'

# Directory containing demultiplexed ONT FASTQs (output of 00_demux_htpamda.py)
ONT_FASTQ_DIR = './demux_results'

# Illumina untreated library controls (already processed via the official
# HT-PAMDA Illumina pipeline; provides the t=0 baseline -- see Methods 1.13)
ILLUMINA_CONTROL_LIB1 = './illumina_control/lib1_raw_counts.csv.gz'
ILLUMINA_CONTROL_LIB2 = './illumina_control/lib2_raw_counts.csv.gz'
ILLUMINA_CONTROL_SAMPLE_LIB1 = 'QC1_LIB1'
ILLUMINA_CONTROL_SAMPLE_LIB2 = 'QC2_LIB2'

# PAM parameters (Methods 1.12, 1.15)
PAM_ORIENTATION = 'three_prime'
MAX_PAM_LENGTH = 8   # full 8xN window sequenced
PAM_START = 4        # positions 0-3 from spacer end (first 4 of 8N)
PAM_LENGTH = 4        # core window reported in heatmaps (Methods 1.15)

# Timepoints in SECONDS: 0=control(t0), 60=1min, 480=8min, 1920=32min (Methods 1.9, 1.14)
TIMEPOINTS = [0, 60, 480, 1920]

# Spacer sequences (Methods 1.9: SPACER1 example given as
# GGGCACGGGCAGCTTGCCGG for HT-PAMDA library 1)
SPACERS = {
    'SPACER1': 'GGGCACGGGCAGCTTGCCGG',
    'SPACER2': 'REPLACE_WITH_LIBRARY_2_SPACER',
}

# Normalization and rate fitting parameters (Methods 1.13, 1.14)
TOP_N_NORMALIZE = 5
INIT_RATE_EST = [0.0001, 0.001, 0.01]
READ_SUM_MIN = 4
TPS_SUM_MIN = 1
USE_TIMEPOINTS = None  # None = use all

# Heatmap parameters (Methods 1.15)
AVERAGE_SPACER = True
HEATMAP_FIXED_MIN = -5.0
HEATMAP_FIXED_MAX = -1.5
LOG_SCALE_HEATMAP = True
PAM1_NT_RANK = {1: 'A', 2: 'C', 3: 'G', 4: 'T'}
PAM2_NT_RANK = {1: 'A', 2: 'C', 3: 'G', 4: 'T'}
PAM1_INDEX_RANK = None
PAM2_INDEX_RANK = None

# ── SAMPLE DEFINITIONS ──────────────────────────────────────────────────
# (sample_id, enzyme_description, library_number, replicate, timepoint_min)
# Replace with the actual candidate/enzyme names for your run. Sample IDs
# must match filenames produced by 00_demux_htpamda.py (SAMPLES list there).
ONT_SAMPLES = [
    ('1A1', 'CANDIDATE_1', 1, 1, 1),
    ('8A1', 'CANDIDATE_1', 1, 1, 8),
    ('32A1', 'CANDIDATE_1', 1, 1, 32),
    # ... add remaining candidates/replicates/timepoints
]

# ══════════════════════════════════════════════════════════════════════════


def main():
    print('HT-PAMDA Custom Pipeline -- ONT + Illumina control')
    print('Run: ' + RUN_NAME)
    print('PAM: orientation={} start={} length={}'.format(
        PAM_ORIENTATION, PAM_START, PAM_LENGTH))
    print('Timepoints (s): {}'.format(TIMEPOINTS))

    barcode_csv = generate_heatmaps.make_barcode_csv(ONT_SAMPLES, RUN_NAME)

    raw_count_files = extract_pam_counts.ont_fastq2count(
        RUN_NAME, ONT_FASTQ_DIR, ONT_SAMPLES, SPACERS,
        TIMEPOINTS, MAX_PAM_LENGTH)

    for lib, control_csv, control_sample in [
        (1, ILLUMINA_CONTROL_LIB1, ILLUMINA_CONTROL_SAMPLE_LIB1),
        (2, ILLUMINA_CONTROL_LIB2, ILLUMINA_CONTROL_SAMPLE_LIB2),
    ]:
        if lib not in raw_count_files:
            print('\nNo data for Library {} -- skipping'.format(lib))
            continue

        lib_spacers = {'SPACER{}'.format(lib): SPACERS['SPACER{}'.format(lib)]}
        run_name_lib = '{}_lib{}'.format(RUN_NAME, lib)

        print('\n' + '=' * 60)
        print('Library {} -- spacer: SPACER{}'.format(lib, lib))
        print('Illumina control: ' + control_sample)

        os.makedirs('output/{}'.format(run_name_lib), exist_ok=True)
        shutil.copy(raw_count_files[lib],
                    'output/{}/PAMDA_1_raw_counts.csv.gz'.format(run_name_lib))

        norm_csv = baseline_correction.rawcount2normcount(
            run_name_lib, raw_count_files[lib], control_csv, control_sample,
            PAM_ORIENTATION, PAM_LENGTH, PAM_START, lib_spacers, TIMEPOINTS,
            MAX_PAM_LENGTH, TOP_N_NORMALIZE)

        rate_csv = kinetic_modeling.normcount2rate(
            run_name_lib, PAM_LENGTH, PAM_START, TIMEPOINTS,
            INIT_RATE_EST, READ_SUM_MIN, TPS_SUM_MIN, USE_TIMEPOINTS,
            input_csv=norm_csv)

        lib_samples = set(
            '{}_Lib{}_Rep{}'.format(e, l, r)
            for _, e, l, r, _ in ONT_SAMPLES if l == lib)
        bc_df = pd.read_csv(barcode_csv)
        bc_df_lib = bc_df[bc_df['sample'].isin(lib_samples)]
        bc_path_lib = 'output/{}/ont_samples_lib{}.csv'.format(RUN_NAME, lib)
        bc_df_lib.to_csv(bc_path_lib, index=False)

        generate_heatmaps.rate2heatmap(
            run_name_lib, bc_path_lib, PAM_LENGTH, PAM_START,
            PAM1_NT_RANK, PAM2_NT_RANK, PAM1_INDEX_RANK, PAM2_INDEX_RANK,
            AVERAGE_SPACER, HEATMAP_FIXED_MIN, HEATMAP_FIXED_MAX,
            LOG_SCALE_HEATMAP, input_csv=rate_csv)

    print('\n' + '=' * 60)
    print('COMPLETE. Results in output/ and figures/ subdirectories.')


if __name__ == '__main__':
    main()
