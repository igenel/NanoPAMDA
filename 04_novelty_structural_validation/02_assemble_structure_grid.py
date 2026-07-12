#!/usr/bin/env python3
"""
Assemble individually rendered PyMOL structure images into a single labeled
grid figure with a shared pLDDT legend (Supplementary Materials figure).

Usage:
    python 02_assemble_structure_grid.py <render_dir> <output_prefix> [n_rows] [n_cols]

Produces <output_prefix>.png and <output_prefix>.pdf.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.patches as mpatches
import os
import sys
import math

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python 02_assemble_structure_grid.py <render_dir> <output_prefix> [n_rows] [n_cols]")
        sys.exit(1)

    render_dir = sys.argv[1]
    output_prefix = sys.argv[2]

    files = sorted([f for f in os.listdir(render_dir) if f.endswith('.png')])
    n = len(files)

    if len(sys.argv) >= 5:
        n_rows, n_cols = int(sys.argv[3]), int(sys.argv[4])
    else:
        n_cols = 4
        n_rows = math.ceil(n / n_cols)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4.5 * n_rows))
    axes = axes.flatten()

    for i, fname in enumerate(files):
        ax = axes[i]
        img = mpimg.imread(os.path.join(render_dir, fname))
        ax.imshow(img)
        ax.axis('off')
        label = fname.replace('.png', '')
        if len(label) > 28:
            label = label[:25] + "..."
        ax.set_title(label, fontsize=8, fontfamily='monospace', pad=2)

    for j in range(n, len(axes)):
        axes[j].axis('off')

    legend_elements = [
        mpatches.Patch(facecolor='#1B4B8C', edgecolor='k', label='pLDDT > 90 (very high)'),
        mpatches.Patch(facecolor='skyblue', edgecolor='k', label='70 < pLDDT \u2264 90 (confident)'),
        mpatches.Patch(facecolor='yellow', edgecolor='k', label='50 < pLDDT \u2264 70 (low)'),
        mpatches.Patch(facecolor='orange', edgecolor='k', label='pLDDT \u2264 50 (very low)'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4, fontsize=10, frameon=False, bbox_to_anchor=(0.5, 0.0))

    fig.suptitle("ESMFold-predicted structures", fontsize=13, fontweight='bold', y=0.995)

    plt.tight_layout(rect=[0, 0.03, 1, 0.98])
    plt.savefig(f"{output_prefix}.png", dpi=200, bbox_inches='tight')
    plt.savefig(f"{output_prefix}.pdf", bbox_inches='tight')
    print(f"Saved {output_prefix}.png / .pdf")
