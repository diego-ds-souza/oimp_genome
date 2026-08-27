#!/usr/bin/env bash
#
# Refine the predicted models against the transcript evidence. `funannotate
# update` runs PASA over the models from step 06, adding UTRs, correcting exon
# boundaries and adding alternative transcripts, then quantifies the result
# with kallisto to drop transcripts with no expression support.
#
# Run from the repository root:
#   conda activate oimp_03_annotation
#   bash 03_annotation/07_funannotate_update.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
OUTDIR="${OUTDIR:-results/03_annotation/funannotate}"
THREADS="${THREADS:-50}"
# -----------------------------------------------------------------------------

command -v funannotate >/dev/null || { echo "ERROR: funannotate not found" >&2; exit 1; }
[[ -d "$OUTDIR" ]] || { echo "ERROR: no such directory: $OUTDIR" >&2; exit 1; }

# 1) Update in place; the directory already holds the training and predict output.
echo "[1/1] funannotate update (${THREADS} cpus)"
funannotate update -i "$OUTDIR" --cpus "$THREADS"

echo "wrote ${OUTDIR}/update_results/"
echo "      ${OUTDIR}/update_misc/kallisto.tsv"
