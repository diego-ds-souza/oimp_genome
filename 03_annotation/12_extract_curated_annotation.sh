#!/usr/bin/env bash
#
# Build the curated annotation from the keep list: subset the GFF3 to the kept
# genes and everything below them, standardise it with AGAT, and pull out the
# matching proteins.
#
# The subset is taken in two passes because child features point at transcript
# IDs rather than gene IDs: the first pass collects the transcripts of kept
# genes, the second emits genes, transcripts and their children together.
#
# Run from the repository root:
#   conda activate oimp_03_annotation
#   bash 03_annotation/12_extract_curated_annotation.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
GFF3="${GFF3:-results/03_annotation/funannotate/annotate_results/Oncideres_impluviata.gff3}"
PROTEINS="${PROTEINS:-results/03_annotation/funannotate/annotate_results/Oncideres_impluviata.proteins.fa}"
OUTDIR="${OUTDIR:-results/03_annotation/curated}"
KEEP="${KEEP:-results/03_annotation/curated/keep_genes.txt}"
NAME="${NAME:-Oncideres_impluviata}"
# -----------------------------------------------------------------------------

for t in awk seqkit; do
  command -v "$t" >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }
done
for f in "$GFF3" "$PROTEINS" "$KEEP"; do
  [[ -s "$f" ]] || { echo "ERROR: missing or empty: $f" >&2; exit 1; }
done

mkdir -p "$OUTDIR"
keep=$(mktemp); tx=$(mktemp)
trap 'rm -f "$keep" "$tx"' EXIT
sed 's/\r$//' "$KEEP" | awk 'NF' | sort -u > "$keep"
echo "kept genes: $(wc -l < "$keep")"

# 1) Collect the transcripts belonging to kept genes.
echo "[1/4] collecting transcripts"
awk -F'\t' -v keep="$keep" '
  function attr(s, k,   n,i,a,p) {
    n = split(s, a, /;/)
    for (i = 1; i <= n; i++) { split(a[i], p, /=/); if (p[1] == k) return p[2] }
    return ""
  }
  BEGIN { while ((getline g < keep) > 0) if (g != "") kept[g] = 1 }
  /^#/ { next }
  $3 == "mRNA" || $3 == "tRNA" || $3 == "transcript" {
    if (attr($9, "Parent") in kept) print attr($9, "ID")
  }' "$GFF3" | sort -u > "$tx"
echo "kept transcripts: $(wc -l < "$tx")"

# 2) Emit the subset: kept genes, their transcripts, and any feature whose
#    Parent is one of those transcripts.
echo "[2/4] subsetting the GFF3"
awk -F'\t' -v keep="$keep" -v txlist="$tx" '
  function attr(s, k,   n,i,a,p) {
    n = split(s, a, /;/)
    for (i = 1; i <= n; i++) { split(a[i], p, /=/); if (p[1] == k) return p[2] }
    return ""
  }
  function parent_kept(s,   p,n,i,P) {
    p = attr(s, "Parent"); if (p == "") return 0
    n = split(p, P, /,/)
    for (i = 1; i <= n; i++) if (P[i] in tx) return 1
    return 0
  }
  BEGIN {
    while ((getline g < keep)  > 0) if (g != "") kept[g] = 1
    while ((getline t < txlist) > 0) if (t != "") tx[t] = 1
  }
  /^#/ { print; next }
  $3 == "gene" { if (attr($9, "ID") in kept) print; next }
  $3 == "mRNA" || $3 == "tRNA" || $3 == "transcript" { if (attr($9, "ID") in tx) print; next }
  { if (parent_kept($9)) print }' "$GFF3" > "${OUTDIR}/${NAME}.clean.gff3"

# 3) Standardise attribute order and feature nesting.
echo "[3/4] AGAT"
if command -v agat_convert_sp_gxf2gxf.pl >/dev/null; then
  agat_convert_sp_gxf2gxf.pl -g "${OUTDIR}/${NAME}.clean.gff3" \
                             -o "${OUTDIR}/${NAME}.clean_AGAT.gff3" >/dev/null
else
  echo "NOTE: AGAT not found; keeping the unstandardised GFF3" >&2
  cp "${OUTDIR}/${NAME}.clean.gff3" "${OUTDIR}/${NAME}.clean_AGAT.gff3"
fi

# 4) Subset the proteins to the kept transcripts.
echo "[4/4] extracting proteins"
seqkit grep -f "$tx" "$PROTEINS" > "${OUTDIR}/${NAME}.proteins.clean.fa"

echo "wrote ${OUTDIR}/${NAME}.clean_AGAT.gff3"
echo "      ${OUTDIR}/${NAME}.proteins.clean.fa  ($(grep -c '^>' "${OUTDIR}/${NAME}.proteins.clean.fa") proteins)"
