#!/usr/bin/env python
"""
Custom demultiplexing of single-step PCR amplicon pools for HT-PAMDA
Nanopore sequencing runs (Methods 1.11).

Reads are computationally binned based on their unique 10-nucleotide forward
and reverse index pairs (embedded during the single-step PCR, Section 1.10).
Barcode matching allows up to 2 mismatches (Hamming distance) and checks
both the read sequence and its reverse complement, since Nanopore
sequencing natively sequences both strands.

Usage:
    python 00_demux_htpamda.py input.fastq output_dir [max_mismatches] [search_window]

SAMPLES below defines the barcode-pair -> sample-ID mapping for one example
run; edit this list (or the F/R barcode sequences, matching
Supplementary Table S3) for each new sequencing run.
"""
from __future__ import print_function
import sys, os, gzip

# Example sample sheet for one run: (sample_id, forward_barcode, reverse_barcode)
# Forward barcodes (F1-F9) and reverse barcodes (R1-R3, one per timepoint: 1/8/32 min)
# See Supplementary Table S3 for the full oligo sequences.
SAMPLES = [
    ("1A1",  "ATCACGTTGC", "CTTATGGCTG"),  # F1, R1
    ("1A2",  "CGTGTAAGCG", "CTTATGGCTG"),  # F2, R1
    ("1A3",  "GACTACTGCA", "CTTATGGCTG"),  # F3, R1
    ("1A4",  "CATGGCATTA", "CTTATGGCTG"),  # F4, R1
    ("1A5",  "AGTCTCGCAG", "CTTATGGCTG"),  # F5, R1
    ("1A6",  "TTGTGATGAC", "CTTATGGCTG"),  # F6, R1
    ("1B1",  "GTCGTAGTCA", "CTTATGGCTG"),  # F7, R1
    ("1B2",  "TCTGTATCTC", "CTTATGGCTG"),  # F8, R1
    ("1B3",  "GCTTAGAGAC", "CTTATGGCTG"),  # F9, R1
    ("8A1",  "ATCACGTTGC", "TGACAGTCGC"),  # F1, R2
    ("8A2",  "CGTGTAAGCG", "TGACAGTCGC"),  # F2, R2
    ("8A3",  "GACTACTGCA", "TGACAGTCGC"),  # F3, R2
    ("8A4",  "CATGGCATTA", "TGACAGTCGC"),  # F4, R2
    ("8A5",  "AGTCTCGCAG", "TGACAGTCGC"),  # F5, R2
    ("8A6",  "TTGTGATGAC", "TGACAGTCGC"),  # F6, R2
    ("8B1",  "GTCGTAGTCA", "TGACAGTCGC"),  # F7, R2
    ("8B2",  "TCTGTATCTC", "TGACAGTCGC"),  # F8, R2
    ("8B3",  "GCTTAGAGAC", "TGACAGTCGC"),  # F9, R2
    ("32A1", "ATCACGTTGC", "CAGTAACTCA"),  # F1, R3
    ("32A2", "CGTGTAAGCG", "CAGTAACTCA"),  # F2, R3
    ("32A3", "GACTACTGCA", "CAGTAACTCA"),  # F3, R3
    ("32A4", "CATGGCATTA", "CAGTAACTCA"),  # F4, R3
    ("32A5", "AGTCTCGCAG", "CAGTAACTCA"),  # F5, R3
    ("32A6", "TTGTGATGAC", "CAGTAACTCA"),  # F6, R3
    ("32B1", "GTCGTAGTCA", "CAGTAACTCA"),  # F7, R3
    ("32B2", "TCTGTATCTC", "CAGTAACTCA"),  # F8, R3
    ("32B3", "GCTTAGAGAC", "CAGTAACTCA"),  # F9, R3
]


def rc(seq):
    t = {"A": "T", "T": "A", "G": "C", "C": "G", "a": "t", "t": "a", "g": "c", "c": "g"}
    return "".join(t.get(b, b) for b in reversed(seq))


def hamming(a, b):
    if len(a) != len(b):
        return 999
    return sum(x != y for x, y in zip(a, b))


def find_bc(seq, bc_set, window=100, max_mm=2):
    search = seq[:window].upper()
    bc_len = 10
    best_bc = None
    best_mm = max_mm + 1
    for bc in bc_set:
        for i in range(len(search) - bc_len + 1):
            mm = hamming(search[i:i + bc_len], bc)
            if mm < best_mm:
                best_mm = mm
                best_bc = bc
    return best_bc, best_mm


def demux(seq, fwd_set, rev_set, window=100, max_mm=2):
    """Try both strand orientations since ONT sequences either strand randomly."""
    seq_rc = rc(seq)
    f1, mm1 = find_bc(seq, fwd_set, window, max_mm)
    r1, mm2 = find_bc(seq_rc, rev_set, window, max_mm)
    if f1 and r1:
        return f1, r1
    r2, mm3 = find_bc(seq, rev_set, window, max_mm)
    f2, mm4 = find_bc(seq_rc, fwd_set, window, max_mm)
    if r2 and f2:
        return f2, r2
    return None, None


def main():
    if len(sys.argv) < 3:
        print("Usage: python 00_demux_htpamda.py input.fastq output_dir [max_mm] [window]")
        sys.exit(1)
    infile = sys.argv[1]
    outdir = sys.argv[2]
    max_mm = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    window = int(sys.argv[4]) if len(sys.argv) > 4 else 100

    fwd_set = set(s[1] for s in SAMPLES)
    rev_set = set(s[2] for s in SAMPLES)
    pair_map = dict(((s[1], s[2]), s[0]) for s in SAMPLES)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    handles = {}
    for name, _, _ in SAMPLES:
        handles[name] = open(os.path.join(outdir, name + ".fastq"), "w")
    handles["unassigned"] = open(os.path.join(outdir, "unassigned.fastq"), "w")

    fh = gzip.open(infile, "rt") if infile.endswith(".gz") else open(infile)
    total = assigned = 0
    counts = dict((n, 0) for n in list(handles.keys()))

    while True:
        hdr = fh.readline().rstrip()
        if not hdr:
            break
        seq = fh.readline().rstrip()
        plus = fh.readline().rstrip()
        qual = fh.readline().rstrip()
        total += 1
        if total % 10000 == 0:
            sys.stderr.write("\r" + str(total) + " reads...")
            sys.stderr.flush()

        sample = "unassigned"
        if len(seq) >= 50:
            fbc, rbc = demux(seq, fwd_set, rev_set, window, max_mm)
            if fbc and rbc:
                s = pair_map.get((fbc, rbc))
                if s:
                    sample = s
                    assigned += 1

        handles[sample].write(hdr + "\n" + seq + "\n" + plus + "\n" + qual + "\n")
        counts[sample] += 1

    fh.close()
    for h in handles.values():
        h.close()

    sys.stderr.write("\n")
    print("Total:     " + str(total))
    print("Assigned:  " + str(assigned) + " (" + str(round(100.0 * assigned / total, 1)) + "%)")
    print("\nReads per sample:")
    for name, _, _ in SAMPLES:
        if counts[name] > 0:
            print("  " + name + ": " + str(counts[name]))


if __name__ == "__main__":
    main()
