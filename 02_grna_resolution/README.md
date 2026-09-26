# 02_grna_resolution — Guide RNA architectural resolution and length filtering (Methods 1.3)

Unlike `01_discovery/`, this stage is a mix of scripted steps and **manual
curation**. Both are documented here so the pipeline is fully reproducible
end-to-end; the manual steps are not placeholders awaiting a script — they
were genuinely performed by hand and are recorded as a protocol instead.

## Scripted steps

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `01_extract_candidate_fastas.py` | 20kb-proximity-filtered hits (`01_discovery/09_filter_20kb_proximity.py` output) + renamed contig `.fna` files | Per-candidate `.fasta` files + `results.tsv` manifest |
| 2 | `02_run_crisprtracrrna.sh` | Per-candidate `.fasta` files | CRISPRtracrRNA per-candidate `.csv` outputs (array/anti-repeat/tracrRNA scores) |
| 3 | `03_select_and_truncate_tracr.py` | CRISPRtracrRNA `.csv` outputs + Cas9 candidate TSV (`01_discovery/10_build_cas_tsv.py`) | `tracr_processed_output.tsv`: one selected, terminator-truncated tracrRNA per locus, paired with its Cas9 protein sequence |

**External dependency:** [CRISPRtracrRNA](https://github.com/BackofenLab/CRISPRtracrRNA)
(BackofenLab) must be installed separately, including its own
`CRISPRidentify` and `CRISPRcasIdentifier` dependencies. Run with
`--model_type II`. See that repository's README for full setup instructions.

## Manual curation steps (no script — protocol below)

Performed on the output of step 3 above (`tracr_processed_output.tsv`):

1. **tracrRNA length filtering (70–150 nt).** Applied manually in spreadsheet
   software to the full candidate pool (this is the filter responsible for
   the 11,694 → 7,168 and 865 → 319 reduction reported in Results). Records
   with a `final_tracr` length outside this window were excluded.

2. **Non-stringent quality checks**, applied manually only to candidates
   surviving clustering (`03_clustering/`) and novelty/catalytic-residue
   screening (`04_novelty_structural_validation/`) — none of these were used
   as strict exclusionary filters:
   - Presence of a GYY repeat-start motif
   - Anti-repeat (AR) tail initiation coordinate
   - Minimum one base match within the first 3 nucleotides of the repeat sequence

3. **2D secondary structure visualization**, via the [NuPACK web server](https://www.nupack.org/)
   (default parameters), performed manually for the surviving cluster
   representatives (not the full candidate pool). Structures were manually
   screened for stable hairpin formation, proper nexus architecture, and a
   distinct nexus junction.

4. **sgRNA chimera assembly**, performed manually only for the final 7
   MAG-derived candidates selected for experimental validation: the
   processed repeat was linked to the truncated tracrRNA via a GAAA
   tetraloop, and the resulting 2D structure was validated a final time via
   NuPACK.

2D NuPACK structures for the 7 MAG-derived candidates selected for
experimental validation are provided in Supplementary Materials.
