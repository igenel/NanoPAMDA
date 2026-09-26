#!/usr/bin/env python3
"""
Generate PAM preference heatmaps from fitted rate constants (Methods 1.15).

Mirrors the official HT-PAMDA pipeline's rate2heatmap logic exactly.
Averages rate constants across spacer targets, log-transforms, and plots on
a fixed color scale (-5.0 to -1.5) for direct visual comparison across
different nuclease variants.
"""
import os
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm.autonotebook import tqdm


def make_barcode_csv(ont_samples, run_name):
    """
    Generate a barcode CSV for rate2heatmap from the ONT sample list.
    Used only for sample descriptions in heatmap titles.
    """
    groups = {}
    for sid, enzyme, lib, rep, time_min in ont_samples:
        key = '{}_Lib{}_Rep{}'.format(enzyme, lib, rep)
        desc = '{}_Lib{}_Rep{}'.format(enzyme, lib, rep)
        groups[key] = desc

    bc_path = 'output/{}/ont_samples.csv'.format(run_name)
    os.makedirs('output/{}'.format(run_name), exist_ok=True)
    with open(bc_path, 'w') as f:
        f.write('sample,description,P5_sample_barcode,P7_sample_barcode\n')
        for i, (key, desc) in enumerate(sorted(groups.items())):
            f.write('{},{},AAAA,TTTT\n'.format(key, desc))
    return bc_path


def rate2heatmap(run_name, barcode_csv, pam_length, pam_start,
                  pam1_nt_rank, pam2_nt_rank,
                  pam1_idx_rank=None, pam2_idx_rank=None,
                  avg_spacer=True, heatmap_fixed_min=-5.0, heatmap_fixed_max=-1.5,
                  log_scale=True, input_csv=None):
    """
    Generate heatmap PDFs from rate constants. Averages across spacer
    targets (Methods 1.15). Each entry in 'Sample' (e.g.
    'CANDIDATE_Lib1_Rep1') retains its own library/replicate identity
    throughout -- this function does not average or otherwise combine
    across replicates.
    """
    print('\n=== STEP 4: rate2heatmap ===')
    plt.switch_backend('agg')

    if input_csv is None:
        input_csv = 'output/{}/PAM_start_{}_length_{}/PAMDA_3_rates.csv'.format(
            run_name, pam_start, pam_length)

    variant_ids = pd.read_csv(barcode_csv)
    variant_name_dict = dict(zip(variant_ids['sample'], variant_ids['description']))

    fig_dir = 'figures/{}/PAM_start_{}_length_{}'.format(run_name, pam_start, pam_length)
    os.makedirs(fig_dir, exist_ok=True)

    df_input = pd.read_csv(input_csv)

    if pam1_idx_rank is None or pam2_idx_rank is None:
        split = int(pam_length // 2)
        pam1_idx_rank = list(range(0, split))[::-1]
        pam2_idx_rank = list(range(split, pam_length))
    else:
        split = len(pam1_idx_rank)

    df_input['PAM_pt1'] = df_input['PAM'].str[:split]
    df_input['PAM_pt2'] = df_input['PAM'].str[split:]
    spacers_in_data = df_input['Spacer'].unique().tolist()
    pam_length_actual = len(df_input['PAM'].iloc[0])

    numbers = ['1', '2', '3', '4']
    pam_space = [''.join(p) for p in itertools.product(numbers, repeat=pam_length_actual)]
    pam1_ids = sorted(set(int(p[:split]) for p in pam_space))
    pam2_ids = sorted(set(int(p[split:]) for p in pam_space))

    indices = []
    for pam in pam1_ids:
        s = str(pam)
        tmp = {i: j for i, j in zip(pam1_idx_rank, list(s.zfill(split)))}
        indices.append(''.join(pam1_nt_rank[int(tmp[i])] for i in sorted(tmp)))

    columns = []
    cols_width = pam_length_actual - split
    for pam in pam2_ids:
        s = str(pam)
        tmp = {i: j for i, j in zip(pam2_idx_rank, list(s.zfill(cols_width)))}
        columns.append(''.join(pam2_nt_rank[int(tmp[i])] for i in sorted(tmp)))

    colors = {'A': '#feb2b1', 'C': '#14c7fe', 'T': '#14f485', 'G': '#f8ffa3'}
    scaling_dict = {
        2: {1: [1, 3.5, 4, 1.07]},
        3: {1: [1, 1.7, 1, 2.15], 2: [1, 2, 4, 3.53]},
        4: {1: [1, 0.7, 0.3, 5.5], 2: [1, 1.7, 1, 4.3], 3: [1, 1, 4, 14.1]},
        5: {1: [1, 0.25, 0.08, 17], 2: [1, 0.6, 0.3, 11.5], 3: [1, 1, 1, 14.15], 4: [1, 0.3, 3, 56.5]},
    }

    pbar = tqdm(desc='heatmaps', total=df_input['Sample'].nunique())
    for variant in df_input['Sample'].unique():
        desc = variant_name_dict.get(variant, variant)
        new_cols = ['{}-{}'.format(col, variant) for col in columns]
        df_out = pd.DataFrame(columns=new_cols, index=indices)

        if avg_spacer:
            for row_idx in df_out.index:
                for col_name in df_out.columns:
                    pam2 = col_name.split('-')[0]
                    sample = '-'.join(col_name.split('-')[1:])
                    rate_avg = 0.0
                    for sp in sorted(spacers_in_data):
                        vals = df_input.loc[
                            (df_input['PAM_pt1'] == row_idx) &
                            (df_input['PAM_pt2'] == pam2) &
                            (df_input['Sample'] == sample) &
                            (df_input['Spacer'] == sp), 'Rate_Constant_k'].tolist()
                        if vals:
                            try:
                                rate_avg += float(vals[0])
                            except Exception:
                                pass
                    rate_avg /= max(len(spacers_in_data), 1)
                    df_out.loc[row_idx, col_name] = np.log10(rate_avg) if log_scale else rate_avg

            df_out = df_out.astype(float)

            hmap_min = heatmap_fixed_min if heatmap_fixed_min else None
            hmap_max = heatmap_fixed_max if heatmap_fixed_max else None

            if 2 <= pam_length_actual <= 5:
                axes = [4 * (pam_length_actual - split), 4 * split]
                fig, ax = plt.subplots(1, figsize=(axes[0], axes[1]))
                sns.heatmap(df_out, ax=ax, vmin=hmap_min, vmax=hmap_max,
                            square=True, cmap='Blues', cbar=True,
                            cbar_kws={'shrink': axes[1] / axes[0] / 2,
                                      'label': 'Log10(rate)' if log_scale else 'rate',
                                      'aspect': 8},
                            linewidth=0.2, linecolor='White',
                            xticklabels=False, yticklabels=False)
                hs = fig.get_size_inches()
                sd = scaling_dict.get(pam_length_actual, {}).get(split, [1, 1, 1, 1])

                x_text = [[columns[n][m] for n in range(len(columns))] for m in range(len(columns[0]))][::-1]
                x_color = [[colors[c] for c in row] for row in x_text]
                xt = ax.table(cellText=x_text, cellColours=x_color, cellLoc='center', loc='top')
                for _, cell in xt.get_celld().items():
                    cell.set_linewidth(0)

                y_text = [[indices[n][m] for m in range(len(indices[0]))] for n in range(len(indices))]
                y_color = [[colors[c] for c in row] for row in y_text]
                yt = ax.table(cellText=y_text, cellColours=y_color,
                               colWidths=[0.06] * len(y_text[0]),
                               cellLoc='center', loc='left')
                for _, cell in yt.get_celld().items():
                    cell.set_linewidth(0)
                xt.set_fontsize(8)
                yt.set_fontsize(8)
                xt.scale(sd[0], sd[1])
                yt.scale(sd[2], hs[1] / sd[3])
                plt.tight_layout(pad=6)
            else:
                fig, ax = plt.subplots()
                plt.title('{} ({})'.format(variant, desc), y=1)
                sns.heatmap(df_out, vmin=hmap_min, vmax=hmap_max,
                            square=True, cmap='Blues', cbar=True,
                            cbar_kws={'shrink': 0.5}, linewidth=0.2,
                            linecolor='White', ax=ax)
                fig.tight_layout(pad=4)
                ax.xaxis.tick_top()
                ax.tick_params(length=0)
                plt.yticks(rotation=0)
                plt.xticks(rotation=90)

            pdf_path = '{}/PAMDA_HEATMAP_{}_{}.pdf'.format(fig_dir, variant, desc)
            csv_path = '{}/PAMDA_HEATMAP_{}_{}.csv'.format(fig_dir, variant, desc)
            plt.savefig(pdf_path)
            df_out.to_csv(csv_path)
            plt.close()
            print('  Saved: ' + pdf_path)

        pbar.update()
    pbar.close()
