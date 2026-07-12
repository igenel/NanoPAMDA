# Cas9 Discovery and Characterization

Computational pipeline and analysis code for the discovery, guide RNA resolution, structural
validation, and PAM characterization of E10 and 31 additional candidate Cas9 orthologs mined
from public metagenome-assembled genome (MAG) catalogs and controlled-access dbGaP raw-read
studies, as described in [paper citation / bioRxiv link once available].

## Repository structure

Each directory corresponds to a Methods subsection in the manuscript, in pipeline order:

| Directory | Manuscript section | Description |
|---|---|---|
| `01_discovery/` | 1.1–1.2 | Full discovery pipeline: header standardization → PILER-CR → Prodigal → hmmsearch + filtering → domtblout cleaning/merging → CDS coordinate mapping → 20kb CRISPR-array proximity filtering → final Cas9 candidate MultiFASTA/TSV database. See `01_discovery/README.md` for the per-script breakdown and `run_pipeline.sh` for the full orchestrated run. |
| `02_grna_resolution/` | 1.3 | Anti-repeat mapping (CRISPRone), repeat:anti-repeat pairing (IntaRNA), terminator detection (Arnold), 2D structure (NuPACK), tracrRNA length filtering, sgRNA assembly |
| `03_clustering/` | 1.4 | CD-HIT redundancy reduction (50% MAG arm / 90% dbGaP arm) |
| `04_novelty_structural_validation/` | 1.5 | BLASTp novelty screening, six-reference catalytic residue MSA, domain architecture check, ESMFold structure prediction, Foldseek TM-score comparison |
| `05_pam_prediction/` | 1.6 | Protein2PAM in silico PAM prediction |
| `06_htpamda_analysis/` | 1.9–1.15 | HT-PAMDA demultiplexing, PAM extraction, baseline correction, kinetic modeling, heatmap generation |
| `supplementary_tables/` | — | Links/pointers to Supplementary Data (hosted separately due to file size; see below) |
| `data/` | — | Source dataset manifest (see Supplementary Table S1) |

## Data availability

Raw source data are not stored in this repository due to size. See:
- **Supplementary Table S1** for the full list of 17 public MAG catalogs (with citations and accession counts) and 2 dbGaP raw-read studies used in this work.
- dbGaP-derived raw sequencing data require controlled-access authorization (accessions: phs000228.v4.p1, phs002232.v1.p1) and are not redistributed here.
- Processed candidate tables (tracrRNA-length-filtered, clustered, and final validated candidate sets) are provided as Supplementary Data files: `[Zenodo/institutional repository DOI — to be added]`.
- Raw Nanopore sequencing data (POD5/FASTQ) generated in this study: `[SRA/ENA accession — to be added]`.

## Requirements

See `environment.yml` for the full computational environment. Key dependencies:
- PILER-CR v1.06
- Prodigal v2.6.3
- HMMER v3.3.2
- IntaRNA v3.3.1
- Arnold (terminator finder)
- NuPACK (default parameters)
- CD-HIT
- Clustal Omega
- ESMFold
- Foldseek
- Protein2PAM
- Dorado v1.4.0 (Oxford Nanopore basecaller)
- Python 3.10+ (Biopython, SciPy, pandas, NumPy)

## Citation

If you use this pipeline or data, please cite: [full citation once published]

## License

Code released under the MIT License (see `LICENSE`). Sequence data and candidate annotations
released under CC-BY 4.0 where applicable — see individual data file headers.

## Contact

[corresponding author email]
