# 04_novelty_structural_validation — Novelty screening, catalytic residue validation, and structural analysis (Methods 1.5)

Most of this stage was performed **manually via web-based tools**, not
scripted. This README documents the protocol actually followed; only the
structure prediction and rendering steps involved code.

## Manual protocol (no scripts — web tools)

Applied to the 424 MAG-derived cluster representatives (Section 1.4) unless
otherwise noted:

1. **Novelty screening (MAG arm only).** Each representative sequence was
   searched against the NCBI non-redundant (nr) protein database using the
   [NCBI web BLAST server](https://blast.ncbi.nlm.nih.gov/) (blastp).
   Candidates sharing ≥90% amino acid identity with any known protein were
   excluded. *(dbGaP arm: not applied.)*

2. **2D gRNA structural filtering (MAG arm only).** Candidates surviving
   novelty screening were required to satisfy the stringent 2D gRNA
   structural rules described in `02_grna_resolution/README.md` (stable
   hairpin, proper nexus architecture, distinct nexus junction), as assessed
   via the NuPACK web server. *(dbGaP arm: not applied.)*

3. **Catalytic residue conservation (both arms).** Surviving candidates
   (MAG arm) or all 67 cluster representatives (dbGaP arm, screened
   sequentially until 16 passing candidates were found — see main Methods)
   were aligned against a six-reference Cas9 panel (SpCas9, FnCas9,
   Nme2Cas9, Cas9, SauCas9, CjCas9) using the
   [Clustal Omega EMBL-EBI web server](https://www.ebi.ac.uk/jdispatcher/msa/clustalo).
   A candidate was retained only if it possessed an allowable residue, using
   SpCas9 numbering, at each of 14 catalytic positions (D10; S15/N15/Y15;
   R66; R70/K70; R74/T74/D74; R78/L78; R165/E165; E762/D762; H840; N854;
   N863; H982/S982; H983; D986) in at least one of the six reference frames.

4. **Domain architecture check (MAG arm only).** The same Clustal Omega
   alignment was visually inspected for large deletions or missing segments
   within known catalytic domains; candidates with such deletions were
   excluded.

5. **Structural similarity via Foldseek (MAG arm only).** Predicted
   structures (see below) were compared against reference Cas9 structures
   using the [Foldseek web server](https://search.foldseek.com/search).
   A TM-score threshold of 0.4 was used to distinguish structurally
   conserved candidates.

## Scripted steps

| Step | Script | Status |
|---|---|---|
| Structure prediction | *(ESMFold script)* | **Pending — to be added** |
| Cartoon rendering + pLDDT coloring | `01_render_structures_pymol.py` | Done |
| Grid figure assembly | `02_assemble_structure_grid.py` | Done |

Once the ESMFold script is added, this stage's script chain will be:
predict structure (ESMFold) → render (`01_render_structures_pymol.py`) →
assemble supplementary grid figure (`02_assemble_structure_grid.py`).

`01`/`02` were used to generate the Supplementary Figure showing all 16
MAG-derived candidates selected for experimental validation, colored by
pLDDT (AlphaFold-style 4-tier scheme). Requires PyMOL (open-source build:
`pip install pymol-open-source`) and matplotlib.

Structures are rendered independently per candidate (not superimposed on a
common reference), since several candidates are highly sequence-divergent
from SpCas9 (~30% identity in some cases) and a common structural alignment
was not attempted.
