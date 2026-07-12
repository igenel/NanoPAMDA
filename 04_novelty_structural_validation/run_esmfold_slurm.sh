#!/bin/bash
# SLURM job script used internally to run 00_run_esmfold.py on an HPC
# cluster via a Singularity/Apptainer container (Methods 1.5).
#
# NOTE: the container image referenced here (esmfold.sif) is a large,
# cluster-specific build and is NOT redistributed with this repository.
# It is not required to reproduce this analysis -- 00_run_esmfold.py can be
# run directly against a standard `esm` installation (see that script's
# docstring and https://github.com/facebookresearch/esm for setup
# instructions). This job script is included only to document exactly how
# structure prediction was executed for this study.
#
#SBATCH -J esmfold
#SBATCH --output=esmfold_log.txt
#SBATCH --error=esmfold_err.txt
#SBATCH -p gpu2dq
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --gres=gpu:1
#SBATCH --mem=128G
#SBATCH -t 01:00:00

# ---------------------------------------------------------------------------
# CONFIG - adjust for your environment
# ---------------------------------------------------------------------------
SIF_CONTAINER="${SIF_CONTAINER:-./esmfold.sif}"   # not provided in this repo
MODEL_CACHE_DIR="${MODEL_CACHE_DIR:-./torch_hub_cache}"
INPUT_FASTA="${INPUT_FASTA:-./candidates.fasta}"
OUTPUT_DIR="${OUTPUT_DIR:-./structures}"
# ---------------------------------------------------------------------------

singularity exec --nv --writable-tmpfs \
    "$SIF_CONTAINER" \
    python 00_run_esmfold.py \
    -m "$MODEL_CACHE_DIR" \
    -i "$INPUT_FASTA" \
    -o "$OUTPUT_DIR" --chunk-size=256
