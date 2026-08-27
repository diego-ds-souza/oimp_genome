#!/usr/bin/env python3
"""
Drop gene models resting on ab initio prediction alone, keeping every model
with transcript or protein evidence (the flags from step 11).

The filter is strict: no ab-initio-only model is rescued, however plausible it
looks. A rescue rule was trialled (minimum protein length, start codon, no
internal stops, multi-exon or UTR) and would have kept 1,755 more genes; it was
not used for the published annotation.

Run from the repository root:
    conda activate oimp_03_annotation
    python 03_annotation/11_filter_abinitio_only.py
"""
import csv
import os

# --------------------------------- settings ----------------------------------
OUTDIR = os.environ.get("OUTDIR", "results/03_annotation/curated")
EVIDENCE = os.path.join(OUTDIR, "gene_evidence.tsv")
KEEP_LIST = os.path.join(OUTDIR, "keep_genes.txt")
SUMMARY = os.path.join(OUTDIR, "filter_summary.tsv")
# -----------------------------------------------------------------------------


def flag(row, name):
    return str(row.get(name, "")).strip().lower() in ("true", "t", "1", "yes", "y")


# 1) Classify each gene.
rows = list(csv.DictReader(open(EVIDENCE), delimiter="\t"))
abinitio_only = [r for r in rows
                 if flag(r, "has_abinitio")
                 and not flag(r, "has_rna")
                 and not flag(r, "has_protein")]
dropped = {r["gene_id"] for r in abinitio_only}
kept = [r["gene_id"] for r in rows if r["gene_id"] not in dropped]

# 2) Write the keep list, which steps 13 onwards subset against.
os.makedirs(OUTDIR, exist_ok=True)
with open(KEEP_LIST, "w") as out:
    out.write("\n".join(kept) + "\n")

# 3) Record what the filter did.
with open(SUMMARY, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["total_genes", "abinitio_only", "dropped", "kept", "rescue_enabled"])
    w.writerow([len(rows), len(abinitio_only), len(dropped), len(kept), False])

print(f"kept {len(kept)}/{len(rows)} genes (dropped {len(dropped)} ab-initio-only)")
print(f"wrote {KEEP_LIST}")
print(f"      {SUMMARY}")
