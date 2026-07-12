# 05_pam_prediction — In silico PAM prediction (Methods 1.6)

No script — performed manually via the Protein2PAM web interface.

## Protocol

1. The full-length amino acid sequence of each of the 32 candidates selected
   for experimental validation (16 MAG-derived + 16 dbGaP-derived; see
   `04_novelty_structural_validation/` and Section 1.5) was submitted
   individually to the [Protein2PAM](https://github.com/Profluent-AI/protein2pam)
   web interface (Profluent Bio; Nayfach, Bhatnagar, Novichkov, et al.).
2. Protein2PAM's pretrained 650M-parameter protein language model (pLM)
   encoder + MLP head returns target position-specific PAM nucleotide
   probabilities for each submitted sequence.
3. Predicted PAM probability outputs for each candidate were downloaded and
   compiled into the Cas9 candidate metadata table (see supplementary
   materials).

## Note

This step was applied only to the 32 final candidates, not to the full
cluster-representative pool (424 MAG + 67 dbGaP) — it is a downstream
prediction used to complement the experimental HT-PAMDA characterization
(Sections 1.9–1.15), not a filtering/screening criterion.
