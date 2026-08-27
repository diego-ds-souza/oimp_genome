#!/usr/bin/env python3
"""
Connect the CAFE result to the annotated host-interface gene families
(Table S10b).

For each family it reports how many orthogroups carrying a member of that
family were tested, how many evolved at a rate inconsistent with the
genome-wide background, how many of those expanded on the O. impluviata
branch specifically, and how many genes that branch gained.

It also asks whether orthogroups containing an annotated family member are
significant more often than orthogroups in general, with a one-sided Fisher
exact test.

Run from the repository root:
    conda activate oimp_05_gene_family_evolution
    python 05_gene_family_evolution/07_summarise_cafe.py
"""
import collections
import csv
import os
from math import exp, lgamma

# --------------------------------- settings ----------------------------------
RES = os.environ.get("RES", "results/05_gene_family_evolution")
MODEL = os.environ.get("MODEL", "base_error")     # the run reported in the paper

FAMILY_RESULTS = os.path.join(RES, MODEL, "Base_family_results.txt")
CHANGE_TAB = os.path.join(RES, MODEL, "Base_change.tab")
BRANCH_P = os.path.join(RES, MODEL, "Base_branch_probabilities.tab")
ORTHOGROUPS = os.environ.get(
    "ORTHOGROUPS", "results/04_orthology/orthofinder_out/Orthogroups.tsv")
ASSIGNMENTS = os.environ.get(
    "ASSIGNMENTS", "results/07_transcriptomics/families/gene_family_assignments.tsv")
OUTDIR = os.environ.get("OUTDIR", os.path.join(RES, "summary"))

FOCAL = os.environ.get("FOCAL", "Oncideres_impluviata")
ALPHA = float(os.environ.get("ALPHA", 0.05))

# The 16 host-interface families reported in the paper, in table order.
GROUPS = [
    ("Chemosensory", ["OBP", "CSP", "OR", "GR", "IR", "SNMP"]),
    ("Detoxification", ["P450", "GST", "UGT", "CCE", "ABC"]),
    ("Digestion", ["SerPro", "CysPro", "AspPro", "MetPro", "CAZyme"]),
]
# -----------------------------------------------------------------------------

TARGET = [f for _, fams in GROUPS for f in fams]
GROUP_OF = {f: g for g, fams in GROUPS for f in fams}


def fisher_greater(a, b, c, d):
    """One-sided Fisher exact test, no SciPy dependency."""
    def lchoose(n, k):
        return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    n, r1, c1 = a + b + c + d, a + b, a + c
    denom = lchoose(n, c1)
    p = sum(exp(lchoose(r1, x) + lchoose(n - r1, c1 - x) - denom)
            for x in range(a, min(r1, c1) + 1))
    odds = (a * d) / (b * c) if b and c else float("inf")
    return odds, min(1.0, p)


def focal_column(header):
    """CAFE tags tip columns as Name<n>; find the focal species' column."""
    for i, h in enumerate(header):
        if h.split("<")[0] == FOCAL:
            return i
    raise SystemExit(f"{FOCAL} not found in {header[:5]}...")


# 1) Family-wide p-values: which orthogroups evolved atypically fast.
pval = {}
with open(FAMILY_RESULTS) as f:
    for row in f:
        if row.startswith("#"):
            continue
        og, p = row.split("\t")[:2]
        pval[og] = float(p)

# 2) Change and branch probability on the focal branch, so an expansion can be
#    attributed to O. impluviata rather than to the tree as a whole.
def read_focal(path):
    out = {}
    with open(path) as f:
        header = next(f).rstrip("\n").split("\t")
        col = focal_column(header)
        for row in f:
            c = row.rstrip("\n").split("\t")
            if len(c) > col:
                out[c[0]] = float(c[col])
    return out

change = read_focal(CHANGE_TAB)
branch_p = read_focal(BRANCH_P)

# 3) Map each orthogroup to the annotated families of its O. impluviata genes.
fam_of_gene = {}
with open(ASSIGNMENTS) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        fams = (r.get("families") or "").strip()
        if fams and fams != "NA":
            fam_of_gene[r["gene_id"]] = [x for x in fams.split(",") if x in TARGET]

og_families = collections.defaultdict(set)
with open(ORTHOGROUPS) as f:
    header = next(f).rstrip("\n").split("\t")
    try:
        col = next(i for i, h in enumerate(header) if h.startswith(FOCAL.split("_")[0][:5])
                   or FOCAL.replace("_", "") in h.replace("_", ""))
    except StopIteration:
        raise SystemExit(f"no {FOCAL} column in {ORTHOGROUPS}")
    for row in f:
        c = row.rstrip("\n").split("\t")
        if len(c) <= col:
            continue
        for prot in c[col].split(", "):
            gene = prot.strip().split("|")[-1].split("-T")[0]
            for fam in fam_of_gene.get(gene, []):
                og_families[c[0]].add(fam)

# 4) One row per family, counting only orthogroups CAFE actually tested.
os.makedirs(OUTDIR, exist_ok=True)
rows = []
for fam in TARGET:
    ogs = [og for og, fams in og_families.items() if fam in fams and og in pval]
    sig = [og for og in ogs if pval[og] < ALPHA]
    exp_ = [og for og in sig
            if change.get(og, 0) > 0 and branch_p.get(og, 1) < ALPHA]
    rows.append({
        "Functional group": GROUP_OF[fam],
        "Gene family": fam,
        "Orthogroups tested": len(ogs),
        "Rapidly evolving": len(sig),
        "Expanded on the O. impluviata branch": len(exp_),
        "Genes gained on that branch": int(sum(change[og] for og in exp_)),
    })

# Totals count DISTINCT orthogroups: one orthogroup carries two families and
# would otherwise be counted twice.
all_ogs = {og for og, fams in og_families.items() if fams and og in pval}
all_sig = {og for og in all_ogs if pval[og] < ALPHA}
all_exp = {og for og in all_sig
           if change.get(og, 0) > 0 and branch_p.get(og, 1) < ALPHA}
rows.append({
    "Functional group": "Total", "Gene family": "Total",
    "Orthogroups tested": len(all_ogs), "Rapidly evolving": len(all_sig),
    "Expanded on the O. impluviata branch": len(all_exp),
    "Genes gained on that branch": int(sum(change[og] for og in all_exp)),
})

out_file = os.path.join(OUTDIR, "cafe_family_summary.tsv")
with open(out_file, "w", newline="") as out:
    w = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

# 5) Are annotated orthogroups significant more often than the rest?
tested = set(pval)
a = len(all_sig)
b = len(all_ogs) - a
c = len({og for og in tested - all_ogs if pval[og] < ALPHA})
d = len(tested - all_ogs) - c
odds, p = fisher_greater(a, b, c, d)

print(f"{len(all_ogs)} orthogroups carry an annotated family member; "
      f"{a} rapidly evolving, {len(all_exp)} expanded on the focal branch, "
      f"{int(sum(change[og] for og in all_exp))} genes gained")
print(f"enrichment vs the other {len(tested) - len(all_ogs)} tested "
      f"orthogroups: odds ratio {odds:.2f}, p = {p:.2g}")
print(f"wrote {out_file}")
