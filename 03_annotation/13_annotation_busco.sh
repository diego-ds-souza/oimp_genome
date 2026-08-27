#!/usr/bin/env bash
#
# Completeness of the annotated protein set, before and after filtering
# (Table S8). Run in protein mode against insecta_odb10, so it measures the
# gene set rather than the assembly.
#
# Run from the repository root:
#   conda activate oimp_03_annotation
#   bash 03_annotation/13_annotation_busco.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
INITIAL="${INITIAL:-results/03_annotation/funannotate/annotate_results/Oncideres_impluviata.proteins.fa}"
CURATED="${CURATED:-results/03_annotation/curated/Oncideres_impluviata.proteins.clean.fa}"
OUTDIR="${OUTDIR:-results/03_annotation/busco}"
LINEAGE="${LINEAGE:-insecta_odb10}"
THREADS="${THREADS:-16}"
# -----------------------------------------------------------------------------

command -v busco >/dev/null || { echo "ERROR: busco not found" >&2; exit 1; }
[[ -f "$INITIAL" ]] || { echo "ERROR: no such file: $INITIAL" >&2; exit 1; }

INITIAL=$(readlink -f "$INITIAL")
[[ -f "$CURATED" ]] && CURATED=$(readlink -f "$CURATED")
mkdir -p "$OUTDIR"; cd "$OUTDIR"

# 1) The full annotated set.
echo "[1/2] BUSCO, initial annotation"
busco --in "$INITIAL" --out busco_initial \
      --lineage_dataset "$LINEAGE" --mode proteins --cpu "$THREADS"

# 2) The curated set. Completeness should hold, showing the filter removed
#    unsupported models rather than real genes.
if [[ -f "$CURATED" ]]; then
  echo "[2/2] BUSCO, curated annotation"
  busco --in "$CURATED" --out busco_curated \
        --lineage_dataset "$LINEAGE" --mode proteins --cpu "$THREADS"
else
  echo "[2/2] skipped: $CURATED not found" >&2
fi

echo "wrote ${OUTDIR}/busco_*/short_summary*.txt"
