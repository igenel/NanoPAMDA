#!/usr/bin/env python3
"""
Fit pseudo-first-order exponential decay to normalized PAM counts to
determine the cleavage rate constant per PAM (Methods 1.14).

y(t) = a * exp(-k*t)

Mirrors the official HT-PAMDA pipeline's normcount2rate logic exactly.
Sequences failing the minimum read/timepoint thresholds are excluded and
assigned a null ('NaN') rate.
"""
import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from tqdm.autonotebook import tqdm


def normcount2rate(run_name, pam_length, pam_start, timepoints,
                    init_rate_est=(0.0001, 0.001, 0.01),
                    read_sum_minimum=4,
                    tps_sum_minimum=1,
                    use_timepoints=None, input_csv=None):
    """
    Fit exponential decay y(t) = a * exp(-k*t) to get rate constant k per PAM.
    """
    print('\n=== STEP 3: normcount2rate ===')

    if input_csv is None:
        input_csv = 'output/{}/PAM_start_{}_length_{}/PAMDA_2_norm_counts.csv'.format(
            run_name, pam_start, pam_length)
    df = pd.read_csv(input_csv)

    if use_timepoints is None:
        use_timepoints = list(range(len(timepoints)))

    def func(x, a, b):
        return a * np.exp(-b * x)

    timepoint_indices = [[i, timepoints[i]] for i in use_timepoints]
    ks = []

    pbar = tqdm(desc='samples', total=df['Sample'].nunique())
    prev_sample = None
    for _, row in df.iterrows():
        if row['Sample'] != prev_sample:
            pbar.update()
            prev_sample = row['Sample']

        tps = [x[1] for x in timepoint_indices]
        obs_raw = [1e-11] + [row.get('Raw_Counts_{}'.format(x[0]), 0)
                              for x in timepoint_indices if str(x[0]) != '0']
        obs_norm = [row.get('Norm_Counts_{}'.format(x[0]), 0) for x in timepoint_indices]

        # Remove zero-read timepoints
        zero_idx = [i for i, e in enumerate(obs_raw) if e == 0]
        for zi in sorted(zero_idx, reverse=True):
            del tps[zi]
            del obs_norm[zi]
            del obs_raw[zi]

        # Minimum read/timepoint thresholds (Methods 1.11: sequences failing
        # a minimum of 4 total reads or missing valid counts in >=1
        # timepoint are excluded and assigned a null value)
        if sum(obs_raw) >= read_sum_minimum and len(obs_norm) >= tps_sum_minimum:
            min_search = []
            for j in init_rate_est:
                try:
                    popt, _ = curve_fit(func, tps, obs_norm, p0=[1.0, j])
                    pred = [func(t, popt[0], popt[1]) for t in tps]
                    error = sum((p - o) ** 2 for p, o in zip(pred, obs_norm))
                    min_search.append([error, list(popt)])
                except Exception:
                    continue
            ks.append(sorted(min_search)[0][1][1] if min_search else 'NaN')
        else:
            ks.append('NaN')
    pbar.close()

    valid_ks = [k for k in ks if isinstance(k, float) and k > 0]
    min_k = min(valid_ks) if valid_ks else 1e-6
    df['Rate_Constant_k'] = [k if (isinstance(k, float) and k > 0) else
                              (min_k if k != 'NaN' else 'NaN') for k in ks]

    out_dir = 'output/{}/PAM_start_{}_length_{}'.format(run_name, pam_start, pam_length)
    out_path = '{}/PAMDA_3_rates.csv'.format(out_dir)
    df.to_csv(out_path, index=False)
    print('  Written: ' + out_path)
    return out_path
