# 01_discovery — CRISPR array and Cas9 candidate discovery (Methods 1.1–1.2)

Scripts in this directory are numbered in execution order. `run_pipeline.sh`
runs steps 01–10 end-to-end against a single input dataset directory; each
step can also be run independently using the paths shown below.

| Step | Script | Input | Output |
|---|---|---|---|
| 0 | `00_metahipmer2_assembly.sh` | Raw reads (dbGaP arm only: HMP-A, TOPDECC) | Assembled contigs |
| 1 | `01_rename_headers.py` | Raw/assembled contig FASTA directory | Renamed `.fna` files (`>contig-100_N` headers) |
| 2 | `02_run_pilercr.sh` | Renamed `.fna` files | PILER-CR output + array coordinate summaries |
| 3 | `03_run_prodigal.sh` | PILER-CR array coordinates + renamed contigs | `.faa` (proteins) / `.gbk` (coordinates) per dataset |
| 4 | `04_run_hmmsearch_filter.sh` | `.faa` files + Cas HMM profile library | Filtered `.domtblout` hits (coverage >70%, bitscore ≥40, e-value <1e-5) |
| 5 | `05_clean_domtblout.py` | Filtered `.domtblout` hits | Cleaned contig/protein ID pairs |
| 6 | `06_merge_domtblout.py` | Cleaned hit files | Deduplicated, merged hit lists |
| 7 | `07_get_protein_locations.py` | Merged hits + `.gbk` files | Protein CDS genomic coordinates |
| 8 | `08_join_contig_protein.sh` | Merged hits + CDS coordinates | Joined contig/protein/location table |
| 9 | `09_filter_20kb_proximity.py` | Joined table + PILER-CR array coordinates | Hits within ±20 kb of a CRISPR array |
| 10 | `10_build_cas_tsv.py` | 20kb-filtered hits + `.faa`/`.gbk` files | Final Cas9 candidate MultiFASTA + TSV database |

## Notes on cleanup from internal lab versions

These scripts were adapted from internal cluster-execution versions for public
release:
- Hardcoded personal paths (e.g., `/okyanus/users/.../miniconda3`, personal
  `pilercr` binary paths) were replaced with `$PATH` lookups or environment
  variables (`PILERCR_BIN`, etc.) — see comments in each script.
- SLURM job directives (account, email, partition) were removed from the
  orchestrator; adapt `run_pipeline.sh` to your own scheduler/environment if
  running on a cluster.
- Conda environment activation calls were removed; activate the environment
  from `environment.yml` (or your own equivalent) before running.

## Requirements

`pilercr`, `prodigal`, and `hmmsearch` must be available on `$PATH` (or
pointed to via the relevant environment variable). See the top-level
`environment.yml`.
