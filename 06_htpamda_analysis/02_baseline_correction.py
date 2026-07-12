#!/usr/bin/env python3
"""
Normalize ONT raw PAM counts using the Illumina untreated library as the
t=0 baseline, with top-N enriched-PAM correction (Methods 1.13).

Mirrors the official HT-PAMDA pipeline's rawcount2normcount logic exactly,
adapted to merge in a separately-sequenced Illumina control library as the
t0 timepoint (since the ONT runs themselves only cover t=1/8/32 min).
"""
import os
import itertools
import warnings
import numpy as np
import pandas as pd
from scipy.stats import linregress


def rawcount2normcount(run_name, raw_count_csv, control_rawcount_csv,
                        control_sample, pam_orientation, pam_length, pam_start,
                        spacers, timepoints, max_pam_length=8, top_n=5):
    """
    Normalize ONT raw counts using the Illumina untreated library as t0.
    """
    print('\n=== STEP 2: rawcount2normcount ===')
    warnings.filterwarnings('ignore', category=RuntimeWarning)

    nucleotides = ['A', 'T', 'C', 'G']
    total_pam_space = [''.join(p) for p in itertools.product(nucleotides, repeat=pam_length)]

    df_input = pd.read_csv(raw_count_csv)

    # Merge with Illumina control (provides t0)
    df_control = pd.read_csv(control_rawcount_csv)
    df_input = pd.concat([df_input, df_control], sort=False)
    control_sample_timepoint_fastq = 1  # control Raw_Counts_1 = t0

    # Truncate 8N PAM to the selected core window (Methods 1.15: 4-nt core)
    print('  grouping counts by PAM bases [{}:{}]'.format(pam_start, pam_start + pam_length))
    column_sort = {'Sample': 'first', 'Spacer': 'first'}
    count_columns = df_input.columns.values[3:]
    new_columns = ['Sample', 'Spacer', 'PAM']
    for col in count_columns:
        column_sort[col] = 'sum'
        new_columns.append(col)

    if pam_orientation == 'three_prime':
        df_input['selected_PAM'] = df_input['PAM'].str[pam_start:pam_start + pam_length]
    else:
        df_input['selected_PAM'] = df_input['PAM'].str[
            max_pam_length - pam_length - pam_start: max_pam_length - pam_start]

    df_list = []
    for sample in df_input['Sample'].unique():
        for spacer in df_input['Spacer'].unique():
            temp = df_input[(df_input['Sample'] == sample) & (df_input['Spacer'] == spacer)] \
                .groupby(['selected_PAM'], as_index=False).agg(column_sort)
            df_list.append(temp)
    df = pd.concat(df_list)
    df = df.rename(columns={'selected_PAM': 'PAM'})
    df = df.loc[:, new_columns].reset_index(drop=True)

    # Normalize each timepoint column by total reads per sample+spacer
    for i in range(1, len(timepoints)):
        col = 'Raw_Counts_{}'.format(i)
        if col in df.columns:
            df['Norm_Counts_{}'.format(i)] = \
                df[col] / df.groupby(['Sample', 'Spacer'])[col].transform('sum')

    # Extract t0 from Illumina control
    print('  setting t0 from Illumina control: ' + control_sample)
    control_dict = {sp: {pam: 0 for pam in total_pam_space} for sp in spacers}
    for _, row in df.iterrows():
        if row['Sample'] == control_sample:
            sp = row['Spacer']
            pam = row['PAM']
            if sp in control_dict and pam in control_dict[sp]:
                control_dict[sp][pam] = row['Norm_Counts_{}'.format(
                    control_sample_timepoint_fastq)]

    df['Norm_Counts_0'] = df.apply(
        lambda row: control_dict.get(row['Spacer'], {}).get(row['PAM'], 0), axis=1)

    # Remove control sample
    df = df[df['Sample'] != control_sample].reset_index(drop=True)

    # Find top-N enriched PAMs per sample (baseline/enrichment correction,
    # Methods 1.13: "top 5 most highly enriched PAM sequences ... using
    # linear regression"; median fractional increase used as the
    # normalization scalar)
    print('  finding top-{} enriched PAMs per sample'.format(top_n))
    uptrends = {}
    x = list(range(len(timepoints)))
    for _, row in df.iterrows():
        key = '{}_{}'.format(row['Sample'], row['Spacer'])
        y = [row.get('Norm_Counts_{}'.format(i), 0) for i in range(len(timepoints))]
        slope = linregress(x, y)[0]
        uptrends.setdefault(key, []).append([slope, y])

    uptrend_corrections = {}
    for key, vals in uptrends.items():
        vals_sorted = sorted(vals)
        top_n_vals = [v[1] for v in vals_sorted[-top_n:]]
        uptrend_corrections[key] = [np.median(col) for col in zip(*top_n_vals)]

    # Apply normalization: correct for enrichment, normalize to t0
    print('  normalizing counts')
    for idx, row in df.iterrows():
        key = '{}_{}'.format(row['Sample'], row['Spacer'])
        corrections = uptrend_corrections.get(key, [1.0] * len(timepoints))
        for i in range(len(timepoints)):
            norm_col = 'Norm_Counts_{}'.format(i)
            if norm_col in df.columns:
                val = row.get(norm_col, 0)
                corr = corrections[i] if corrections[i] != 0 else 1e-10
                t0 = row.get('Norm_Counts_0', 0)
                t0 = t0 if t0 != 0 else 1e-10
                df.at[idx, norm_col] = (val / corr) / t0

    out_dir = 'output/{}/PAM_start_{}_length_{}'.format(run_name, pam_start, pam_length)
    os.makedirs(out_dir, exist_ok=True)
    out_path = '{}/PAMDA_2_norm_counts.csv'.format(out_dir)
    df.to_csv(out_path, index=False)
    print('  Written: ' + out_path)
    return out_path
