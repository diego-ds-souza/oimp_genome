#!/usr/bin/env bash
#
# Synteny between O. impluviata and Monochamus alternatus (GCA_037114965.1),
# used to corroborate the chromosome-scale assignment of scaffolds (Fig. 1b).
#
# Alignments are made against the unmasked O. impluviata assembly, keeping hits
# longer than 100 bp with MAPQ >= 30. That leaves roughly 70,000 alignments,
# far too many to plot, so nearby links are bundled: two or more links within a
# 1 Mb window become one connection, which gives the 771 links in the figure.
#
# A coarser bundling (500 or more links within 10 Mb, giving 24 connections)
# was also run for a broader view; set MAX_GAP and MIN_BUNDLE to reproduce it.
#
# Run from the repository root:
#   conda activate oimp_02_assembly
#   bash 02_assembly/10_synteny.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
REFERENCE="${REFERENCE:-results/02_assembly/fixed/Hifiasm_purged.sorted.fixed.renumbered.fasta}"
OUTGROUP="${OUTGROUP:-data/GCA_037114965.1_Malternatus.fna}"   # M. alternatus
OUTDIR="${OUTDIR:-results/02_assembly/synteny}"
THREADS="${THREADS:-8}"

MIN_LEN="${MIN_LEN:-100}"          # minimum alignment length (bp)
MIN_MAPQ="${MIN_MAPQ:-30}"         # minimum mapping quality
MAX_GAP="${MAX_GAP:-1000000}"      # bundlelinks: merge links within 1 Mb
MIN_BUNDLE="${MIN_BUNDLE:-2}"      # bundlelinks: minimum links per bundle
# -----------------------------------------------------------------------------

command -v minimap2 >/dev/null || { echo "ERROR: minimap2 not found" >&2; exit 1; }
for f in "$REFERENCE" "$OUTGROUP"; do
  [[ -f "$f" ]] || { echo "ERROR: no such file: $f" >&2; exit 1; }
done

mkdir -p "$OUTDIR"

# 1) Align M. alternatus against the unmasked O. impluviata assembly.
#    -x asm5 assumes roughly 5% sequence divergence.
echo "[1/3] minimap2"
minimap2 -x asm5 -t "$THREADS" "$REFERENCE" "$OUTGROUP" > "${OUTDIR}/Malt_vs_Oimp.paf"

# 2) Filter, then write Circos link format:
#    outgroup chrom/start/end  <tab>  O. impluviata scaffold/start/end.
#    The MAPQ filter keeps only alignments with a unique best placement.
echo "[2/3] filtering (>= ${MIN_LEN} bp, MAPQ >= ${MIN_MAPQ})"
awk -v L="$MIN_LEN" -v Q="$MIN_MAPQ" 'BEGIN{OFS="\t"}
  $11 >= L && $12 >= Q { print $1,$3,$4,$6,$8,$9 }' \
  "${OUTDIR}/Malt_vs_Oimp.paf" > "${OUTDIR}/links.txt"

# 3) Bundle nearby links so the figure shows syntenic blocks rather than every
#    individual hit.
echo "[3/3] bundlelinks"
command -v bundlelinks >/dev/null || {
  echo "ERROR: bundlelinks not found. It ships with Circos as" >&2
  echo "       tools/bundlelinks/bin/bundlelinks; add that directory to PATH." >&2
  exit 1
}
bundlelinks -links "${OUTDIR}/links.txt" \
            -max_gap "$MAX_GAP" \
            -min_bundle_membership "$MIN_BUNDLE" \
            > "${OUTDIR}/links_bundled.txt"

echo "wrote ${OUTDIR}/links_bundled.txt"
