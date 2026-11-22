# 07_bacterial_assay — E10 ATc-inducible killing assay in *E. coli*

Corresponds to the Methods section "Bacterial Nuclease Activity Assay" and to manuscript
Fig. 4 / Supplementary Table S15.

E10 was expressed from an ATc-inducible promoter (pdCas9-bacteria backbone, Addgene #44249,
dCas9 replaced by E10) together with sgRNAs (psgRNA, Addgene #114005, E10-compatible scaffold)
targeting the essential genes *rpoB* (rpoB-1, rpoB-4) or *gyrA* (gyrA-1, gyrA-2). Double-strand
breaks at essential chromosomal genes are lethal in *E. coli*, so cleavage is read out as loss
of colony formation with ATc. Controls: E10 without a guide RNA and a guide RNA without E10.

| File | Description |
|---|---|
| `cfu_counts.csv` | Raw colony counts for three independent transformations per condition, with dilution and plated volume |
| `plot_fig4.py` | Converts counts to CFU/mL, applies Welch's t-test on log10 CFU, and draws Fig. 4 |

```bash
python plot_fig4.py cfu_counts.csv Figure4_E10_ATc_killing
```

CFU/mL = colonies / (dilution x plated volume). Plates with no colonies are set to the limit of
detection (25 CFU/mL) for targeting guides, so log10 reductions are lower bounds.

An earlier experiment with guides targeting *adhE*, *pta* and *xylB* (n = 2) is reported in
Supplementary Fig. 3 and Supplementary Table S16 of the manuscript.
