#!/usr/bin/env python3
"""
Pull the class-by-divergence matrix out of the RepeatMasker .divsum file.

A .divsum has two parts: a per-family Kimura table, and a matrix giving the
bases covered by each repeat class at each divergence bin. The repeat landscape
plot needs only the second, as a plain space-separated table.

Run from the repository root:
    conda activate oimp_03_annotation
    python 03_annotation/03_divsum_table.py
"""
import os
import sys

# --------------------------------- settings ----------------------------------
OUTDIR = os.environ.get("OUTDIR", "results/03_annotation/repeats")
NAME = os.environ.get("NAME", "Oimp")

DIVSUM = os.path.join(OUTDIR, f"{NAME}.divsum")
OUT_FILE = os.path.join(OUTDIR, f"{NAME}.divsum.table")

MARKER = "Coverage for each repeat class and divergence"
# -----------------------------------------------------------------------------

# 1) Find the matrix header and copy everything after it.
with open(DIVSUM) as f:
    lines = f.read().splitlines()

start = next((i for i, l in enumerate(lines) if MARKER in l), None)
if start is None:
    sys.exit(f"ERROR: '{MARKER}' not found in {DIVSUM}")

matrix = [l for l in lines[start + 1:] if l.strip()]

# 2) Write it out. The first row is the header (Div, then one repeat class per
#    column); every later row is one divergence bin.
with open(OUT_FILE, "w") as out:
    out.write("\n".join(matrix) + "\n")

print(f"wrote {OUT_FILE}  ({len(matrix) - 1} divergence bins, "
      f"{len(matrix[0].split()) - 1} repeat classes)")
