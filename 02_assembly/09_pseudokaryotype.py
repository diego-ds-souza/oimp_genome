#!/usr/bin/env python3
"""
Pseudo-karyotype of the chromosome-scale scaffolds, with telomeric ends marked
(base layer of Fig. 1a; the published figure was finished in Illustrator).

Run from the repository root:
    conda activate oimp_02_assembly
    python 02_assembly/09_pseudokaryotype.py
"""
import csv
import os

import matplotlib.pyplot as plt

# --------------------------------- settings ----------------------------------
OUTDIR = os.environ.get("OUTDIR", "results/02_assembly/tidk")
SUMMARY_FILE = os.path.join(OUTDIR, "telomere_summary.tsv")
OUT_PREFIX = os.path.join(OUTDIR, "Oimp_pseudokaryotype_telomeres")
MIN_LENGTH = int(os.environ.get("MIN_LENGTH", 5_000_000))   # plot scaffolds >= 5 Mb
# -----------------------------------------------------------------------------

# 1) Read the telomere calls, keeping only the chromosome-scale scaffolds.
scaffolds = []
with open(SUMMARY_FILE) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        length = int(row["length_bp"])
        if length < MIN_LENGTH:
            continue
        scaffolds.append({
            "name": row["scaffold"],
            "length": length,
            "t5": row["telomere_5prime"].lower() == "yes",
            "t3": row["telomere_3prime"].lower() == "yes",
        })

scaffolds.sort(key=lambda s: s["length"], reverse=True)

# 2) One bar per scaffold, drawn to scale, with a dot at each telomeric end.
fig, ax = plt.subplots(figsize=(6, 0.4 * len(scaffolds) + 2))
for y, sc in enumerate(scaffolds):
    length_mb = sc["length"] / 1e6
    ax.hlines(y, 0, length_mb, linewidth=4, color="black")
    if sc["t5"]:
        ax.plot(0, y, marker="o", markersize=6, color="firebrick")
    if sc["t3"]:
        ax.plot(length_mb, y, marker="o", markersize=6, color="firebrick")
    ax.text(-0.5, y, sc["name"], va="center", ha="right", fontsize=8)

# 3) Save.
ax.set_xlabel("Scaffold length (Mb)")
ax.set_yticks([])
ax.invert_yaxis()
ax.set_title("Pseudo-karyotype of $\\it{O.\\ impluviata}$ scaffolds with telomeres")

plt.tight_layout()
plt.savefig(f"{OUT_PREFIX}.pdf")
plt.savefig(f"{OUT_PREFIX}.png", dpi=300)
print(f"wrote {OUT_PREFIX}.pdf / .png")
