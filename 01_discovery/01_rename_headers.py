#!/usr/bin/env python3
"""
Standardize FASTA contig headers across source assemblies (Methods 1.1).

Renames every contig header in each .fna/.fasta/.fa file in a directory to
the form '>contig-100_N' (N = sequential index within the file), and
normalizes filenames (underscores/dots -> hyphens) for consistent downstream
parsing across the 17 MAG catalogs and 2 dbGaP-derived assemblies.

This is an idempotent version: files already in the correct header/filename
format are skipped. Files that ARE renamed have their original preserved
alongside as 'old_<original_filename>' rather than being deleted, so no data
is lost if the script is re-run or something goes wrong downstream.

Usage:
    python 01_rename_headers.py <directory_path>
"""
import os
import sys


def change_headers_and_filenames(directory):
    for filename in os.listdir(directory):
        if filename.endswith(('.fna', '.fasta', '.fa')):
            filepath = os.path.join(directory, filename)

            with open(filepath, 'r') as file:
                lines = file.readlines()

            # Check if headers and filenames are already in the correct format
            all_correct = True
            contig_counter = 1
            for line in lines:
                if line.startswith('>'):
                    expected_header = f">contig-100_{contig_counter}\n"
                    if line != expected_header:
                        all_correct = False
                        break
                    contig_counter += 1

            base_name, ext = os.path.splitext(filename)
            expected_filename = base_name.replace('_', '-').replace('.', '-') + ".fna"
            if expected_filename != filename:
                all_correct = False

            if all_correct:
                print(f"{filename} is already in the correct format. Skipping...")
                continue

            # Modify headers and filenames if not in the correct format
            new_lines = []
            contig_counter = 1
            for line in lines:
                if line.startswith('>'):
                    new_header = f">contig-100_{contig_counter}\n"
                    new_lines.append(new_header)
                    contig_counter += 1
                else:
                    new_lines.append(line)

            new_base_name = base_name.replace('_', '-').replace('.', '-')
            new_filename = f"{new_base_name}.fna"
            new_filepath = os.path.join(directory, new_filename)

            with open(new_filepath, 'w') as file:
                file.writelines(new_lines)

            old_filepath = os.path.join(directory, f"old_{filename}")
            os.rename(filepath, old_filepath)
            print(f"Renamed {filename} to {new_filename} and original file to old_{filename}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python 01_rename_headers.py <directory_path>")
        sys.exit(1)
    directory_path = sys.argv[1]
    change_headers_and_filenames(directory_path)
