#!/usr/bin/env python3
"""
Render ESMFold-predicted structures as cartoon representations colored by
per-residue pLDDT confidence (Methods 1.5, structural visualization).

Uses PyMOL (open-source build; `pip install pymol-open-source`) to render
each PDB file individually with AlphaFold-style 4-tier pLDDT coloring
(B-factor field), auto-oriented per structure.

Usage:
    python 01_render_structures_pymol.py <pdb_dir> <output_png_dir>

Note: structures are NOT superimposed on a common reference -- each is
independently auto-oriented by PyMOL for best individual framing. Given the
high sequence divergence among candidates (down to ~30% identity to SpCas9
in some cases), a common structural alignment was not attempted.
"""
import pymol
from pymol import cmd
import os
import sys

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python 01_render_structures_pymol.py <pdb_dir> <output_png_dir>")
        sys.exit(1)

    pdb_dir = sys.argv[1]
    out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    pymol.finish_launching(['pymol', '-qc'])

    pdb_files = sorted([f for f in os.listdir(pdb_dir) if f.endswith('.pdb')])
    print(f"Found {len(pdb_files)} PDB files")

    for pdb_file in pdb_files:
        name = pdb_file.replace('.pdb', '')
        cmd.reinitialize()
        cmd.load(os.path.join(pdb_dir, pdb_file), "struct")
        cmd.hide("everything")
        cmd.show("cartoon", "struct")
        cmd.bg_color("white")

        # AlphaFold-style pLDDT coloring (4-tier, based on B-factor column)
        cmd.color("orange", "struct")             # pLDDT <= 50 (very low)
        cmd.color("yellow", "struct and b > 50")   # 50 < pLDDT <= 70 (low)
        cmd.color("skyblue", "struct and b > 70")  # 70 < pLDDT <= 90 (confident)
        cmd.color("marine", "struct and b > 90")   # pLDDT > 90 (very high)

        cmd.orient("struct")
        cmd.zoom("struct", buffer=3)
        cmd.set("ray_opaque_background", 0)
        cmd.set("cartoon_fancy_helices", 1)
        cmd.png(os.path.join(out_dir, f"{name}.png"), width=900, height=900, dpi=150, ray=1)
        print(f"Rendered {name}")

    print("ALL DONE")
