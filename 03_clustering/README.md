# 03_clustering — Sequence clustering and redundancy reduction (Methods 1.4)

| Step | Script | Input | Output |
|---|---|---|---|
| 0 | `00_assign_unique_ids.py` | Merged, length-filtered candidate TSV (MAG-derived or dbGaP-derived) | Same TSV + sequential `unique_id` column (`seq_1`, `seq_2`, ...) |
| 1 | `01_build_multifasta.py` | `unique_id`-assigned candidate TSV | MultiFASTA with `>seq_N` headers |
| 2a | `02_run_cdhit_mag_50pct.sh` | MAG-arm multiFASTA | CD-HIT representatives (50% identity) + `.clstr` file |
| 2b | `02_run_cdhit_dbgap_90pct.sh` | dbGaP-arm multiFASTA | CD-HIT representatives (90% identity) + `.clstr` file |
| 3 | `03_extract_cluster_map.sh` | `.clstr` file | Two-column `(unique_id, cluster_id)` map |
| 4 | `04_merge_cluster_assignments.py` | `unique_id`-assigned candidate TSV + cluster map | Candidate TSV with cluster assignment column added |

Representative sequence per cluster: **CD-HIT's default** (the longest
sequence in the cluster), not a randomly selected one — Methods 1.4 has been
corrected to reflect this.

`unique_id` (`seq_1`, `seq_2`, ...) is assigned once, after all candidate
data for an arm has been collected and length-filtered, and is independent
of the `natural_<i>_<cas_class>` style labels used elsewhere for
readability/tracking. It is the identifier CD-HIT clusters on and that the
cluster-map extraction/merge steps key on throughout this stage.

