#!/usr/bin/env bash
#
# InterProScan over the updated protein set. Run separately rather than through
# `funannotate iprscan` so the XML can be reused; step 09 reads it directly.
#
# Run from the repository root:
#   conda activate oimp_03_annotation
#   bash 03_annotation/08_interproscan.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
PROTEINS="${PROTEINS:-results/03_annotation/funannotate/update_results/Oncideres_impluviata.proteins.fa}"
OUTDIR="${OUTDIR:-results/03_annotation/interproscan}"
THREADS="${THREADS:-50}"
# -----------------------------------------------------------------------------

command -v interproscan.sh >/dev/null || { echo "ERROR: interproscan.sh not found" >&2; exit 1; }
[[ -f "$PROTEINS" ]] || { echo "ERROR: no such file: $PROTEINS" >&2; exit 1; }

mkdir -p "$OUTDIR"
name=$(basename "$PROTEINS")

# 1) Scan. -goterms and -pa add the GO and pathway cross-references that
#    funannotate carries into the annotation; -dp disables the lookup service
#    so the run is reproducible against the local databases.
echo "[1/1] interproscan (${THREADS} cpus)"
interproscan.sh -i "$PROTEINS" \
                -d "$OUTDIR" \
                -f XML \
                -goterms -pa -dp \
                -cpu "$THREADS"

echo "wrote ${OUTDIR}/${name}.xml"
