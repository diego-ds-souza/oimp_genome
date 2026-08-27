#!/usr/bin/env bash
#
# Extract the gene-to-GO mapping from the eggNOG-mapper annotations. This is
# the background against which every study set in step 00 is tested: all
# O. impluviata genes carrying at least one GO term, with identifiers
# normalised to match the study sets.
#
# Run from the repository root:
#   conda activate oimp_06_go_enrichment
#   bash 06_go_enrichment/01_gene2go_from_eggnog.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
EMAP_PATH="${EMAP_PATH:-data/eggnog.emapper.annotations}"   # .tsv or .tsv.gz
OUT_DIR="${OUT_DIR:-results/06_go_enrichment}"              # output directory
# -----------------------------------------------------------------------------

mkdir -p "$OUT_DIR"

echo "============================================================================="
echo "Gene-to-GO mapping extraction"
echo "============================================================================="
echo "Input:  ${EMAP_PATH}"
echo "Output: ${OUT_DIR}/Oi_gene2go.tsv"
echo ""

python3 - << 'PY' "$EMAP_PATH" "$OUT_DIR"
import sys, os, re, gzip

emap, out_dir = sys.argv[1:3]
out = os.path.join(out_dir, "Oi_gene2go.tsv")

def opf(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p, "r")

GO_COL_CANDS = ("GOs","GO")   # typical eggNOG headers

def norm_id(gid: str) -> str:
    """
    Normalize gene IDs to match study sets:
    - drop species prefix if present: 'O_impluviata_protein|FUN_000001-T1' -> 'FUN_000001-T1'
    - drop isoform suffix: 'FUN_000001-T1' -> 'FUN_000001'
    """
    core = gid.split("|")[-1]
    core = re.sub(r'-T\d+$','', core)         # remove -T1, -T2, ...
    core = re.sub(r'\.t\d+$','', core, flags=re.I) # also handle .t1 if present
    return core

pairs = set()
with opf(emap) as f:
    header = None
    go_idx = None
    for ln in f:
        if ln.startswith("##"):       # meta lines
            continue
        if ln.startswith("#query") or ln.startswith("query"):
            header = ln.lstrip("#").rstrip("\n").split("\t")
            # find GO column by name
            for i, c in enumerate(header):
                if c in GO_COL_CANDS:
                    go_idx = i
                    break
            continue
        if header is None or go_idx is None:
            # haven't seen the real header yet
            continue
        cols = ln.rstrip("\n").split("\t")
        if len(cols) <= go_idx:
            continue
        gid = cols[0]
        gos = cols[go_idx]
        if not gos or gos == "-":
            continue
        gene = norm_id(gid)
        for go in re.split(r"[;, ]+", gos):
            go = go.strip()
            if go.startswith("GO:"):
                pairs.add((gene, go))

with open(out, "w") as g:
    for gene, go in sorted(pairs):
        g.write(f"{gene}\t{go}\n")

genes = len({p[0] for p in pairs})
go_terms = len({p[1] for p in pairs})
print(f"Genes with >=1 GO annotation: {genes}")
print(f"Unique GO terms: {go_terms}")
print(f"Total gene-GO pairs: {len(pairs)}")
print(f"\nWrote: {out}")
PY

echo ""
echo "DONE: ${OUT_DIR}/Oi_gene2go.tsv"
