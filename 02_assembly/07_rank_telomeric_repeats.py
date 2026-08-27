#!/usr/bin/env python3
"""
Rank the curated Coleoptera telomeric motifs by how often they occur, using the
per-window counts from `tidk find` in 06_telomere_repeats.sh.

This is what identified AACCT (TTAGG) as the motif to carry forward: it was by
far the most abundant of the four queried, ahead of ACCTG and AACCC.

Run from the repository root:
    conda activate oimp_02_assembly
    python 02_assembly/07_rank_telomeric_repeats.py
"""
import csv
import os
from collections import defaultdict

# --------------------------------- settings ----------------------------------
OUTDIR = os.environ.get("OUTDIR", "results/02_assembly/tidk")
PREFIX = os.environ.get("PREFIX", "Oimp")

WINDOWS = os.path.join(OUTDIR, f"{PREFIX}_find_telomeric_repeat_windows.tsv")
OUT_FILE = os.path.join(OUTDIR, "frequent_telomeric_repeats.csv")
# -----------------------------------------------------------------------------

# 1) Sum forward and reverse counts across every window, per motif.
totals = defaultdict(int)
with open(WINDOWS) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        totals[row["telomeric_repeat"]] += (int(row["forward_repeat_number"])
                                            + int(row["reverse_repeat_number"]))

# 2) Write the motifs in descending order of abundance.
ranked = sorted(totals.items(), key=lambda kv: kv[1], reverse=True)
with open(OUT_FILE, "w", newline="") as out:
    w = csv.writer(out)
    w.writerow(["telomeric_repeat", "total_count"])
    w.writerows(ranked)

for motif, count in ranked:
    print(f"{motif}\t{count}")
print(f"wrote {OUT_FILE}")
