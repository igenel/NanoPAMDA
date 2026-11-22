#!/usr/bin/env python3
"""Plot the E10 ATc killing assay (manuscript Fig. 4) from cfu_counts.csv.

CFU/mL = colonies / (dilution x plated volume). Plates with no colonies are set to the
limit of detection (one colony per plated volume = 25 CFU/mL) for plotting and statistics
(targeting guides only). Welch's t-test is applied to log10-transformed CFU/mL
(No ATc vs +ATc).

Usage: python plot_fig4.py cfu_counts.csv Figure4_E10_ATc_killing
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
from scipy import stats

PINK, BLUE = "#FFA6A4", "#007AF0"
csv, out = sys.argv[1], sys.argv[2]
df = pd.read_csv(csv)
REPS = ["rep1_colonies", "rep2_colonies", "rep3_colonies"]
V = float(df["plated_volume_mL"].iloc[0])
LOD = 1.0 / V


def cfu(row):
    counts = np.array([row[r] for r in REPS], float)
    values = counts / (float(row["dilution_plated"]) * V)
    if row["target_gene"] != "control":
        values = np.where(counts == 0, LOD, values)
    return values


order = list(dict.fromkeys(df["construct"]))
data = {}
for g in order:
    sub = df[df.construct == g]
    a = cfu(sub[sub.condition == "No ATc"].iloc[0])
    b = cfu(sub[sub.condition == "+ATc"].iloc[0])
    data[g] = (a, b)


def stars(p):
    return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < .05 else "ns"


def label(g):
    for gene in ("rpoB", "gyrA"):
        if gene in g:
            return g.replace(gene, r"$\it{%s}$" % gene)
    return g


plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.linewidth": 0.8, "pdf.fonttype": 42})
x = np.arange(len(order)); w = 0.36
fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7.2, 7.4), gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.42})
rng = np.random.default_rng(1)
for i, g in enumerate(order):
    a, b = data[g]
    for off, v, c in ((-w / 2 - 0.01, a, PINK), (w / 2 + 0.01, b, BLUE)):
        m, s = v.mean(), v.std(ddof=1)
        ax.bar(x[i] + off, m, w, color=c, edgecolor="black", lw=0.8, zorder=2)
        ax.errorbar(x[i] + off, m, yerr=[[min(s, m * 0.999)], [s]], color="black", capsize=3, lw=0.9, zorder=3)
        ax.scatter(x[i] + off + rng.uniform(-.07, .07, len(v)), v, s=16, facecolor="white", edgecolor="black", lw=0.8, zorder=4)
    p = stats.ttest_ind(np.log10(a), np.log10(b), equal_var=False).pvalue
    top = max(a.max(), b.max()) * 2.5
    ax.plot([x[i] - w / 2, x[i] - w / 2, x[i] + w / 2 + .01, x[i] + w / 2 + .01], [top * .8, top, top, top * .8], color="black", lw=0.8)
    ax.text(x[i], top * 1.15, stars(p), ha="center", va="bottom", fontsize=9)
ax.set_yscale("log"); ax.set_ylim(10, 3e8)
ax.set_yticks([10 ** k for k in range(1, 9)]); ax.yaxis.set_major_formatter(matplotlib.ticker.LogFormatterSciNotation())
ax.set_ylabel("Colony count (CFU/mL)"); ax.set_xticks(x); ax.set_xticklabels([label(g) for g in order], rotation=30, ha="right")
ax.spines[["top", "right"]].set_visible(False)
ax.bar(0, 0, color=PINK, edgecolor="black", label="No ATc"); ax.bar(0, 0, color=BLUE, edgecolor="black", label="+ATc")
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0, 1.02))
ax.text(-0.12, 1.06, "A", transform=ax.transAxes, fontsize=14, fontweight="bold")
for i, g in enumerate(order):
    a, b = data[g]
    r = np.log10(a) - np.log10(b); m, s = r.mean(), r.std(ddof=1)
    ax2.bar(x[i], m, 0.55, color=BLUE, edgecolor="black", lw=0.8, zorder=2)
    ax2.errorbar(x[i], m, yerr=s, color="black", capsize=3, lw=0.9, zorder=3)
    ax2.scatter(x[i] + rng.uniform(-.1, .1, len(r)), r, s=16, facecolor="white", edgecolor="black", lw=0.8, zorder=4)
    ax2.text(x[i], max(m + s, r.max()) + 0.25, ("≥" if "control" != df[df.construct == g].target_gene.iloc[0] else "") + f"{abs(m) if abs(m) < 0.05 else m:.1f}", ha="center", va="bottom", fontsize=8)
ax2.set_ylim(-0.5, 6.6); ax2.axhline(0, color="black", lw=0.8)
ax2.set_ylabel("log$_{10}$ reduction\n(No ATc / +ATc)"); ax2.set_xticks(x); ax2.set_xticklabels([label(g) for g in order], rotation=30, ha="right")
ax2.spines[["top", "right"]].set_visible(False)
ax2.text(-0.12, 1.06, "B", transform=ax2.transAxes, fontsize=14, fontweight="bold")
for ext, kw in (("png", {"dpi": 600}), ("pdf", {})):
    fig.savefig(f"{out}.{ext}", bbox_inches="tight", facecolor="white", **kw)
