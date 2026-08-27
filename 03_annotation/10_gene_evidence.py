#!/usr/bin/env python3
"""
Score every gene model in the annotated GFF3 by the evidence supporting it.

Three flags are recorded per gene:

  has_protein   the transcript carries an InterPro or Pfam cross-reference, an
                eggNOG, COG or MEROPS note, or a GO term - that is, homology or
                domain evidence;
  has_rna       the gene has 5' or 3' UTR features, which funannotate adds only
                where PASA matched a transcript, or names an RNA-derived source
                in its product or note;
  has_abinitio  the model exists at all, which every prediction does.

A gene with only the third flag rests on ab initio prediction alone. Step 12
removes those.

Run from the repository root:
    conda activate oimp_03_annotation
    python 03_annotation/10_gene_evidence.py
"""
import csv
import os
from collections import defaultdict

# --------------------------------- settings ----------------------------------
GFF3 = os.environ.get(
    "GFF3",
    "results/03_annotation/funannotate/annotate_results/Oncideres_impluviata.gff3")
OUTDIR = os.environ.get("OUTDIR", "results/03_annotation/curated")
OUT_FILE = os.path.join(OUTDIR, "gene_evidence.tsv")

PROTEIN_DBXREF = ("InterPro", "PFAM")
PROTEIN_NOTE = ("EggNog", "COG", "MEROPS")
RNA_KEYWORDS = ("rna", "transcript", "pasa", "stringtie", "trinity",
                "rnaseq", "est", "cdna", "expressed", "mrna", "trna")
# -----------------------------------------------------------------------------


def attrs(field):
    out = {}
    for part in field.rstrip(";").split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
    return out


# 1) Walk the GFF3 once, collecting evidence per gene.
genes = {}
order = []
utr = defaultdict(bool)
exons = defaultdict(int)

with open(GFF3) as f:
    for line in f:
        if line.startswith("#"):
            continue
        col = line.rstrip("\n").split("\t")
        if len(col) < 9:
            continue
        kind, a = col[2], attrs(col[8])

        if kind == "gene":
            gid = a.get("ID", "")
            order.append(gid)
            genes[gid] = {"gene_id": gid, "chromosome": col[0],
                          "start": col[3], "end": col[4], "strand": col[6],
                          "has_rna": False, "has_protein": False,
                          "has_abinitio": False, "has_utrs": False,
                          "total_exons": 0, "transcripts": 0}

        elif kind == "mRNA":
            gid = a.get("Parent", "")
            g = genes.get(gid)
            if g is None:
                continue
            g["transcripts"] += 1
            g["has_abinitio"] = True          # every model is a prediction

            dbxref = a.get("Dbxref", "")
            note = a.get("note", "")
            ontology = a.get("Ontology_term", "")
            product = a.get("product", "")

            if any(k in dbxref for k in PROTEIN_DBXREF) \
                    or any(k in note for k in PROTEIN_NOTE) \
                    or "GO:" in ontology:
                g["has_protein"] = True

            blob = (product + " " + note).lower()
            if any(k in blob for k in RNA_KEYWORDS):
                g["has_rna"] = True

        elif kind in ("five_prime_UTR", "three_prime_UTR"):
            utr[a.get("Parent", "").split("-")[0]] = True

        elif kind == "exon":
            exons[a.get("Parent", "").split("-")[0]] += 1

# 2) UTRs are transcript evidence: funannotate adds them from the PASA
#    alignments, so a model with UTRs was matched by RNA-seq.
for gid, g in genes.items():
    if utr.get(gid):
        g["has_rna"] = True
        g["has_utrs"] = True
    g["total_exons"] = exons.get(gid, 0)

# 3) Write one row per gene.
os.makedirs(OUTDIR, exist_ok=True)
cols = ["gene_id", "chromosome", "start", "end", "strand", "transcripts",
        "total_exons", "has_rna", "has_protein", "has_abinitio", "has_utrs"]
with open(OUT_FILE, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for gid in order:
        w.writerow(genes[gid])

only = sum(1 for g in genes.values()
           if g["has_abinitio"] and not g["has_rna"] and not g["has_protein"])
print(f"wrote {OUT_FILE}  ({len(genes)} genes, {only} ab-initio-only)")
