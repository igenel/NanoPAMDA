#!/usr/bin/env python3
"""
Extract PAM counts from demultiplexed ONT FASTQs (Methods 1.12).

For each read, searches for an exact match to the target spacer sequence
(both forward and reverse-complement orientations) and extracts the
immediately adjacent 8-nucleotide PAM window. Reads are grouped into
(enzyme, library, replicate) sample groups, with one count column per
ONT timepoint.

This is a component of the custom ONT-adapted HT-PAMDA pipeline: each step
here mirrors the corresponding step of the official HT-PAMDA pipeline
(Walton et al., 2021; https://github.com/kleinstiverlab/HT-PAMDA) but is
reimplemented to handle single-end Nanopore long reads rather than
paired-end Illumina reads.

See run_pipeline.py for the full orchestrated run and CONFIG section
(spacer sequences, timepoints, sample sheet) that must be edited per run.
"""
import os
import gzip
import itertools
from tqdm.autonotebook import tqdm


def reverse_complement(seq):
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(comp.get(c, 'N') for c in seq.upper()[::-1])


def find_spacer_and_extract_pam(seq, spacer, max_pam_len):
    """
    Find spacer in seq (exact match) and extract the PAM immediately 3'.
    Also tries the reverse complement of seq. Returns a PAM string of
    length max_pam_len, or None if no valid match is found.
    """
    for s in [seq.upper(), reverse_complement(seq)]:
        loc = s.find(spacer)
        if loc != -1:
            pam_start = loc + len(spacer)
            pam_end = pam_start + max_pam_len
            if pam_end <= len(s):
                pam = s[pam_start:pam_end]
                if all(b in 'ACGT' for b in pam):
                    return pam
    return None


def ont_fastq2count(run_name, ont_fastq_dir, ont_samples, spacers,
                     timepoints, max_pam_len=8):
    """
    Read demultiplexed ONT FASTQs (single-end, one file per sample).
    Build raw count tables in the same format as the official pipeline's
    fastq2count output: Sample, Spacer, PAM, Raw_Counts_1, Raw_Counts_2, ...

    ont_samples: list of (sample_id, enzyme, library_number, replicate, timepoint_min)
    spacers: dict of {spacer_name: spacer_sequence}
    """
    print('\n=== STEP 1: ONT fastq2count ===')

    nucleotides = ['A', 'T', 'C', 'G']
    total_pam_space = [''.join(p) for p in itertools.product(nucleotides, repeat=max_pam_len)]
    n_timepoints = len(timepoints) - 1  # exclude t0 (control)

    groups = {}
    for sid, enzyme, lib, rep, time_min in ont_samples:
        group_key = '{}_Lib{}_Rep{}'.format(enzyme, lib, rep)
        spacer_name = 'SPACER{}'.format(lib)
        if group_key not in groups:
            groups[group_key] = {'spacer': spacer_name, 'lib': lib, 'entries': []}
        groups[group_key]['entries'].append((time_min, sid))

    for g in groups:
        groups[g]['entries'].sort()

    output_by_lib = {1: {}, 2: {}}

    pbar = tqdm(desc='groups', total=len(groups))
    for group_key, info in sorted(groups.items()):
        spacer_name = info['spacer']
        spacer_seq = spacers[spacer_name]
        lib = info['lib']
        entries = info['entries']

        pam_counts = {pam: [0] * len(entries) for pam in total_pam_space}

        for tp_idx, (time_min, sid) in enumerate(entries):
            fastq_path = os.path.join(ont_fastq_dir, sid + '.fastq')
            if not os.path.exists(fastq_path):
                print('  WARNING: not found: ' + fastq_path)
                continue

            total = 0
            extracted = 0
            with open(fastq_path) as fh:
                while True:
                    hdr = fh.readline()
                    if not hdr:
                        break
                    seq = fh.readline().rstrip()
                    fh.readline()
                    fh.readline()
                    total += 1
                    pam = find_spacer_and_extract_pam(seq, spacer_seq, max_pam_len)
                    if pam and pam in pam_counts:
                        pam_counts[pam][tp_idx] += 1
                        extracted += 1

            pct = round(100. * extracted / total, 1) if total > 0 else 0
            tqdm.write('  {} t={}min {}: {}/{} reads ({:.1f}%)'.format(
                group_key, time_min, sid, extracted, total, pct))

        output_by_lib[lib][group_key] = {spacer_name: pam_counts}
        pbar.update()
    pbar.close()

    os.makedirs('output/{}'.format(run_name), exist_ok=True)
    output_files = {}

    for lib in [1, 2]:
        if not output_by_lib[lib]:
            continue
        out_path = 'output/{}/PAMDA_1_raw_counts_lib{}.csv.gz'.format(run_name, lib)
        header = ['Sample', 'Spacer', 'PAM'] + \
                 ['Raw_Counts_{}'.format(i + 1) for i in range(n_timepoints)]
        with gzip.open(out_path, 'wt') as f:
            f.write(','.join(header) + '\n')
            for group_key, spacer_dict in sorted(output_by_lib[lib].items()):
                for spacer_name, pam_dict in spacer_dict.items():
                    for pam, counts in pam_dict.items():
                        row = [group_key, spacer_name, pam] + counts
                        f.write(','.join(map(str, row)) + '\n')
        output_files[lib] = out_path
        print('Written: ' + out_path)

    return output_files
