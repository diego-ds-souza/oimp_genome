#!/usr/bin/env python3
"""
Call a telomere present or absent at each scaffold end from the TIDK window
counts written by 06_telomere_repeats.sh.

A telomere is called when a window within TERMINAL_BP of a scaffold end carries
at least REPEAT_THRESHOLD copies of the motif (forward + reverse). Scaffolds
with a telomere at both ends are the chromosome-scale candidates reported in
the paper.

Run from the repository root:
    conda activate oimp_02_assembly
    python 02_assembly/08_telomere_summary.py
"""
import csv
import os

# --------------------------------- settings ----------------------------------
OUTDIR = os.environ.get("OUTDIR", "results/02_assembly/tidk")
PREFIX = os.environ.get("PREFIX", "Oimp")
MOTIF = os.environ.get("MOTIF", "AACCT")

LEN_FILE = os.path.join(OUTDIR, "scaffold_lengths.tsv")
WINDOWS = os.path.join(OUTDIR, f"{PREFIX}_search_{MOTIF}_telomeric_repeat_windows.tsv")
OUT_FILE = os.path.join(OUTDIR, "telomere_summary.tsv")

TERMINAL_BP = int(os.environ.get("TERMINAL_BP", 50_000))       # end region size
REPEAT_THRESHOLD = int(os.environ.get("REPEAT_THRESHOLD", 30))  # min repeats
# -----------------------------------------------------------------------------

# 1) Scaffold lengths, so a window can be placed relative to each end.
with open(LEN_FILE) as f:
    lengths = {}
    for line in f:
        if line.strip():
            name, length = line.split("\t")[:2]
            lengths[name] = int(length)

summary = {
    name: {"length_bp": length, "telomere_5prime": False, "telomere_3prime": False}
    for name, length in lengths.items()
}

# 2) Walk the TIDK windows and flag the terminal ones that clear the threshold.
with open(WINDOWS) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        scaffold = row["id"]
        if scaffold not in summary:
            continue
        if int(row["forward_repeat_number"]) + int(row["reverse_repeat_number"]) \
                < REPEAT_THRESHOLD:
            continue
        pos = int(row["window"])
        if pos <= TERMINAL_BP:
            summary[scaffold]["telomere_5prime"] = True
        if pos >= summary[scaffold]["length_bp"] - TERMINAL_BP:
            summary[scaffold]["telomere_3prime"] = True

# 3) Write the summary, longest scaffold first.
with open(OUT_FILE, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["scaffold", "length_bp", "telomere_5prime", "telomere_3prime"])
    for scaffold, v in sorted(summary.items(),
                              key=lambda x: x[1]["length_bp"], reverse=True):
        w.writerow([scaffold, v["length_bp"],
                    "yes" if v["telomere_5prime"] else "no",
                    "yes" if v["telomere_3prime"] else "no"])

both = sum(1 for v in summary.values()
           if v["telomere_5prime"] and v["telomere_3prime"])
print(f"wrote {OUT_FILE}  ({both} scaffolds with telomeres at both ends)")
