#!/usr/bin/env bash
#
# Sort the purged contigs by length (longest first) and rename them
# scaffold_1 .. scaffold_N. This is the "contigs were sorted by size and
# renamed sequentially" step in the Methods.
#
# Run from the repository root:
#   conda activate oimp_02_assembly
#   bash 02_assembly/04_sort_and_rename.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
ASM="${ASM:-results/02_assembly/purge_dups/Hifiasm_purged.fa}"  # purged assembly
OUTDIR="${OUTDIR:-results/02_assembly}"                         # output directory
OUT="${OUT:-Hifiasm_purged.sorted.fasta}"                       # output filename
# -----------------------------------------------------------------------------

for t in seqkit samtools awk; do
  command -v "$t" >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }
done
[[ -f "$ASM" ]] || { echo "ERROR: no such file: $ASM" >&2; exit 1; }

mkdir -p "$OUTDIR"
tmp="${OUTDIR}/sorted.tmp.fa"

# 1) Sort by sequence length, longest first.
echo "[1/3] sorting by length"
seqkit sort -l -r "$ASM" > "$tmp"

# 2) Rename the headers in that order.
echo "[2/3] renaming scaffold_1 .. scaffold_N"
awk '/^>/ { print ">scaffold_" ++n; next } { print }' "$tmp" > "${OUTDIR}/${OUT}"
rm -f "$tmp"

# 3) Index, and write the length table used by the telomere and synteny steps.
echo "[3/3] indexing"
samtools faidx "${OUTDIR}/${OUT}"
awk 'BEGIN{OFS="\t"} {print $1,$2}' "${OUTDIR}/${OUT}.fai" \
  > "${OUTDIR}/scaffold_lengths_sorted.tsv"

echo "wrote ${OUTDIR}/${OUT}"
echo "      ${OUTDIR}/scaffold_lengths_sorted.tsv"
