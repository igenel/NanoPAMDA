# NanoPAMDA

Computational pipeline and analysis code for the discovery, guide RNA resolution, structural
validation, and PAM characterization of Cas9 orthologs mined from public environmental metagenome
catalogs and controlled-access human microbiome 
studies, as described in Genel, Suren & Kurt, *NanoPAMDA enables rapid PAM characterization of
metagenome-derived Cas9 effectors*, bioRxiv (2026), https://doi.org/10.64898/2026.11.04.764188.

This pipeline combines scripted, automated steps with manual/web-based curation steps at
several stages. Both are documented equally throughout this repository — manual steps are not
placeholders, they are genuinely part of the methodology and are recorded as explicit protocols
in the relevant subdirectory's README.

## Repository structure

Each directory corresponds to a Methods subsection in the manuscript, in pipeline order:

| Directory | Manuscript section | Status | Description |
|---|---|---|---|
| `01_discovery/` | 1.1–1.2 | Complete | Header standardization → PILER-CR → Prodigal → hmmsearch + filtering → domtblout cleaning/merging → CDS coordinate mapping → 20kb CRISPR-array proximity filtering → final Cas9 candidate MultiFASTA/TSV database. Fully scripted; see `01_discovery/README.md` for the per-script breakdown and `run_pipeline.sh` for the full orchestrated run. |
| `02_grna_resolution/` | 1.3 | Complete | tracrRNA extraction and identification via CRISPRtracrRNA (BackofenLab), representative-element selection, terminator-truncation heuristic (scripted); tracrRNA length filtering, GYY/AR-tail motif checks, NuPACK 2D structure visualization, and sgRNA chimera assembly (manual — see `02_grna_resolution/README.md` for the full protocol). |
| `03_clustering/` | 1.4 | Complete | Unique ID assignment → multiFASTA build → CD-HIT clustering (50% environmental arm / 90% human microbiome arm) → cluster map extraction → merge back into candidate metadata. Fully scripted. |
| `04_novelty_structural_validation/` | 1.5 | Complete | Novelty screening (NCBI web BLAST), catalytic residue MSA (EBI Clustal Omega web server), domain architecture check (visual inspection), and Foldseek TM-score comparison (Foldseek web server) — all manual; see `04_novelty_structural_validation/README.md`. Structure prediction via ESMFold (`00_run_esmfold.py`), PyMOL rendering, and supplementary figure grid assembly scripts included and tested. |
| `05_pam_prediction/` | 1.6 | Complete | In silico PAM prediction via the Protein2PAM web interface (Profluent Bio; Nayfach et al., *Nat. Biotechnol.* 2026), applied manually to the 14 final candidates. No script — protocol documented in `05_pam_prediction/README.md`. |
| `06_htpamda_analysis/` | 1.9–1.15 | Complete | Dorado basecalling (HPC/SLURM, sup model) → custom demultiplexing (10-nt barcode pairs, Hamming distance ≤2) → PAM extraction → reference-library-corrected normalization → pseudo-first-order kinetic rate fitting → PAM preference heatmap generation. A custom reimplementation of the official HT-PAMDA pipeline (Walton et al., 2021) adapted for single-end ONT long reads. Fully scripted; see `06_htpamda_analysis/README.md`. |
| `07_bacterial_assay/` | Methods: Bacterial Nuclease Activity Assay | Complete | E10 ATc-inducible killing assay in *E. coli* (Fig. 4): raw colony counts and the script that computes CFU/mL, statistics and the figure. See `07_bacterial_assay/README.md`. |
| `supplementary_tables/` | — | — | Pointers to Supplementary Data (hosted externally due to file size; see below). |
| `data/` | — | — | Source dataset manifest (see Supplementary Table S1). |

## Data availability

Raw source data are not stored in this repository due to size. See:
- **Supplementary Table S1** for the full list of 17 public metagenome-assembled genome catalogs (with citations and genome counts) and 2 dbGaP raw-read studies used in this work.
- dbGaP-derived raw sequencing data require controlled-access authorization (accessions: phs000228.v4.p1, phs002232.v1.p1) and are not redistributed here.
- Processed candidate tables (tracrRNA-length-filtered, clustered, and final candidate sets for the environmental and human microbiome arms) are provided in Supplementary Tables S3–S8 of the preprint (see `supplementary_tables/README.md`).
- Raw Nanopore sequencing data (POD5/FASTQ) generated in this study: NCBI SRA BioProject PRJNA1504502.
- Plasmid maps for the effector expression vectors (including E10) and the modified BPK1520 gRNA vector are provided with the preprint's supplementary files; ONT barcoded oligos and gBlock/gRNA sequences for all 14 candidates are in Supplementary Tables S13 and S12.

## Requirements

See `environment.yml` for the Python computational environment. Several steps rely on external
command-line tools or web servers rather than Python packages:

**Command-line tools** (must be installed/available on `$PATH` separately):
- PILER-CR v1.06
- Prodigal v2.6.3
- HMMER v3.3.2 (`hmmsearch`)
- MetaHipMer2 (dbGaP raw-read assembly only)
- CRISPRtracrRNA ([BackofenLab/CRISPRtracrRNA](https://github.com/BackofenLab/CRISPRtracrRNA)), including its CRISPRidentify and CRISPRcasIdentifier dependencies
- CD-HIT
- Dorado v1.4.0 (Oxford Nanopore basecaller; GPU required)
- ESMFold (`fair-esm` with ESMFold support; CUDA GPU required — see [facebookresearch/esm](https://github.com/facebookresearch/esm))
- PyMOL (open-source build: `pip install pymol-open-source`) — structure rendering only

**Web servers used manually** (no local installation; see the relevant subdirectory README for the exact protocol followed at each):
- NCBI web BLAST (blastp, nr database) — novelty screening
- EMBL-EBI Clustal Omega — catalytic residue MSA
- NuPACK — 2D gRNA secondary structure prediction
- Foldseek web server — structural similarity / TM-score comparison
- Protein2PAM ([Profluent-AI/protein2pam](https://github.com/Profluent-AI/protein2pam)) — PAM prediction

**Python** 3.10+: Biopython, SciPy, pandas, NumPy, matplotlib, seaborn, tqdm, openpyxl

## Citation

If you use this pipeline or data, please cite: Genel I, Suren C, Kurt IC. NanoPAMDA enables rapid PAM characterization of
metagenome-derived Cas9 effectors. bioRxiv 2026. https://doi.org/10.64898/2026.11.04.764188

## License

Code released under the MIT License (see `LICENSE`). Sequence data and candidate annotations
released under CC-BY 4.0 where applicable — see individual data file headers.

## Contact

Ibrahim C. Kurt (ibrahim.kurt@bogazici.edu.tr)
