#!/usr/bin/env python3
"""
Deduplicate and merge cleaned domtblout hit lists across HMM profile batches
(Methods 1.2), producing one merged hit list per dataset.

Usage:
    python 06_merge_domtblout.py <base_dir_from_05_clean_domtblout>

Writes merged files to <base_dir>/merged/{subdirectory}_{name1}_merged.out
"""
import os
import sys

base_dir = sys.argv[1]
merged_dir = os.path.join(base_dir, "merged")

if not os.path.exists(merged_dir):
    os.makedirs(merged_dir)

file_contents = {}

for root, dirs, files in os.walk(base_dir):
    for filename in files:
        if filename.endswith('_filtered.domtblout'):
            subdirectory = os.path.basename(root)
            filepath = os.path.join(root, filename)

            name_parts = filename.split('_')
            if len(name_parts) >= 2:
                name1 = name_parts[-2]
                with open(filepath, 'r') as file:
                    content = file.read()

                key = (subdirectory, name1)
                if key in file_contents:
                    new_lines = set(content.splitlines())
                    file_contents[key].update(new_lines)
                else:
                    file_contents[key] = set(content.splitlines())

for (subdirectory, name1), contents in file_contents.items():
    merged_content = '\n'.join(contents)
    merged_filename = os.path.join(merged_dir, "{}_{}_merged.out".format(subdirectory, name1))
    with open(merged_filename, 'w') as merged_file:
        merged_file.write(merged_content)
