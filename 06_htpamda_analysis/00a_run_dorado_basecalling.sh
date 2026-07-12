#!/bin/bash
# Dorado basecalling of raw Nanopore signal (POD5) using the super-accurate
# (sup) model (Methods 1.11).
#
# Requires the Dorado basecaller (v1.4.0, Oxford Nanopore Technologies) and
# a GPU-accelerated environment. Adjust SLURM directives / paths below for
# your own cluster and job scheduler.
#
# Usage: adjust CONFIG paths below, then submit via your scheduler (e.g.
# `sbatch 00a_run_dorado_basecalling.sh`) or run directly on a GPU node.

set -euo pipefail

# ---------------------------------------------------------------------------
# CONFIG - adjust for your environment
# ---------------------------------------------------------------------------
DORADO_BIN="${DORADO_BIN:-./dorado-1.4.0-linux-x64/bin/dorado}"
POD5_INPUT="${POD5_INPUT:-./pod5/merged.pod5}"
FASTQ_OUTPUT="${FASTQ_OUTPUT:-./basecalled_sup.fastq}"
# ---------------------------------------------------------------------------

echo "Checking input POD5 exists..."
ls -lh "$POD5_INPUT"

echo "Checking dorado binary exists..."
ls -lh "$DORADO_BIN"

echo "Running Dorado basecaller (sup model)..."
"$DORADO_BIN" basecaller \
    --emit-fastq \
    sup \
    "$POD5_INPUT" \
    > "$FASTQ_OUTPUT"

echo "Basecalling complete: $FASTQ_OUTPUT"
