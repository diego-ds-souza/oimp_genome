#!/usr/bin/env bash
#
# De novo repeat library for O. impluviata with RepeatModeler, then combined
# with the curated Coleoptera repeats from Dfam and Repbase into the custom
# library used for masking.
#
# RepeatModeler runs RECON, RepeatScout and LTR_retriever internally and
# classifies the resulting families against Dfam; it takes about a day on
# 40 threads for this genome.
#
# Run from the repository root:
#   conda activate oimp_03_annotation
#   bash 03_annotation/01_repeat_library.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
GENOME="${GENOME:-results/02_assembly/fixed/Hifiasm_purged.sorted.fixed.renumbered.fasta}"
OUTDIR="${OUTDIR:-results/03_annotation/repeats}"
NAME="${NAME:-Oimp}"                    # RepeatModeler database name
CURATED="${CURATED:-data/Coleoptera_Dfam38_Repbase.fa}"   # curated library
THREADS="${THREADS:-40}"
# -----------------------------------------------------------------------------

for t in BuildDatabase RepeatModeler seqkit; do
  command -v "$t" >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }
done
[[ -f "$GENOME" ]] || { echo "ERROR: no such file: $GENOME" >&2; exit 1; }

GENOME=$(readlink -f "$GENOME")
[[ -f "$CURATED" ]] && CURATED=$(readlink -f "$CURATED")
mkdir -p "$OUTDIR"; cd "$OUTDIR"

# 1) Index the genome for RepeatModeler.
echo "[1/3] BuildDatabase"
BuildDatabase -name "$NAME" "$GENOME"

# 2) De novo family discovery and classification. -LTRStruct adds the
#    structural LTR search, and -genomeSampleSizeMax is set to the whole genome
#    so no sampling is applied. Families land in ${NAME}-families.fa.
gsize=$(seqkit stats "$GENOME" | awk 'NR==2 {gsub(/,/,"",$5); print $5}')
echo "[2/3] RepeatModeler (${THREADS} threads, sample size ${gsize})"
RepeatModeler -database "$NAME" -threads "$THREADS" \
              -genomeSampleSizeMax "$gsize" -LTRStruct > run.out

# 3) Combine the de novo families with the curated Coleoptera library.
#    This concatenated file is what RepeatMasker is given in step 02.
echo "[3/3] building the custom library"
if [[ -f "$CURATED" ]]; then
  cat "${NAME}-families.fa" "$CURATED" > custom_db.fa
else
  echo "NOTE: $CURATED not found; using the de novo families alone" >&2
  cp "${NAME}-families.fa" custom_db.fa
fi

echo "wrote ${OUTDIR}/custom_db.fa"
