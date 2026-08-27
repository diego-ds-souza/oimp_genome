#!/usr/bin/env bash
#
# Functional annotation. `funannotate annotate` combines the InterProScan XML
# from step 08 with its own searches against Pfam, eggNOG/COG, UniProtKB,
# MEROPS, dbCAN/CAZy, BUSCO and the Gene Ontology, and writes the annotated
# GFF3, GenBank and protein files.
#
# MEROPS and dbCAN run here as steps inside the workflow, which is where the
# protease and CAZyme assignments used later come from.
#
# Run from the repository root:
#   conda activate oimp_03_annotation
#   bash 03_annotation/09_funannotate_annotate.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
OUTDIR="${OUTDIR:-results/03_annotation/funannotate}"
IPRSCAN="${IPRSCAN:-results/03_annotation/interproscan/Oncideres_impluviata.proteins.fa.xml}"
BUSCO_DB="${BUSCO_DB:-insecta}"
THREADS="${THREADS:-50}"
# -----------------------------------------------------------------------------

command -v funannotate >/dev/null || { echo "ERROR: funannotate not found" >&2; exit 1; }
[[ -d "$OUTDIR"  ]] || { echo "ERROR: no such directory: $OUTDIR" >&2; exit 1; }
[[ -f "$IPRSCAN" ]] || { echo "ERROR: no such file: $IPRSCAN" >&2; exit 1; }

# 1) Annotate.
echo "[1/1] funannotate annotate (${THREADS} cpus)"
funannotate annotate \
  -i "$OUTDIR" \
  --cpus "$THREADS" \
  --busco_db "$BUSCO_DB" \
  --iprscan "$IPRSCAN"

echo "wrote ${OUTDIR}/annotate_results/Oncideres_impluviata.gff3"
echo "      ${OUTDIR}/annotate_results/Oncideres_impluviata.annotations.txt"
