# 06_htpamda_analysis — HT-PAMDA demultiplexing, PAM extraction, kinetics, and heatmaps (Methods 1.11–1.15)

| Step | Script | Methods section | Input | Output |
|---|---|---|---|---|
| 0a | `00a_run_dorado_basecalling.sh` | 1.11 | Raw POD5 signal files | Basecalled FASTQ (sup model) |
| 0b | `00_demux_htpamda.py` | 1.11 | Basecalled ONT FASTQ | Per-sample demultiplexed FASTQs |
| 1 | `01_extract_pam_counts.py` | 1.12 | Demultiplexed FASTQs | Raw PAM count tables (per library) |
| 2 | `02_baseline_correction.py` | 1.13 | Raw counts + Illumina t0 control | Normalized/baseline-corrected counts |
| 3 | `03_kinetic_modeling.py` | 1.14 | Normalized counts | Fitted rate constants (k) per PAM |
| 4 | `04_generate_heatmaps.py` | 1.15 | Rate constants | PAM preference heatmap PDFs + CSVs |

`run_pipeline.py` orchestrates steps 1–4 end-to-end for both spacer
libraries. Run `00a_run_dorado_basecalling.sh` and `00_demux_htpamda.py`
separately first (once per sequencing run), then point `run_pipeline.py`'s
`ONT_FASTQ_DIR` at the resulting demultiplexed output directory.

## Relationship to the official HT-PAMDA pipeline

This is a **custom reimplementation**, not a direct call to the official
HT-PAMDA codebase (Walton et al., 2021). Steps 1–4 each mirror the logic of
the corresponding official pipeline function (`fastq2count`,
`rawcount2normcount`, `normcount2rate`, `rate2heatmap`) exactly, but are
rewritten to operate on **single-end Nanopore long reads** rather than
paired-end Illumina reads, and to merge in a **separately-generated Illumina
untreated control library** as the t=0 baseline (since the ONT runs
themselves only cover the t=1/8/32 min timepoints — see Section 1.9).

## Basecalling (Methods 1.11)

Dorado basecalling (v1.4.0, super-accurate model) was performed on a
GPU-accelerated high-performance computing (HPC) cluster via SLURM
(`00a_run_dorado_basecalling.sh`), not in a Colab notebook — an earlier
draft of this README incorrectly stated otherwise; that has been corrected.

## Configuration

`00_demux_htpamda.py`'s `SAMPLES` list and `run_pipeline.py`'s `CONFIG`
section (spacer sequences, timepoints, Illumina control paths, and
`ONT_SAMPLES` sample sheet) must be edited for each sequencing run — these
are shown with placeholder/example values only. All hardcoded
cluster-specific paths from the original internal version have been
removed.

## Requirements

`numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`, `tqdm`. See top-level
`environment.yml`.
