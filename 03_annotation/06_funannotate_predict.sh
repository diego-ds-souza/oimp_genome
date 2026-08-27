#!/usr/bin/env bash
#
# Consensus gene models. `funannotate predict` trains and runs Augustus,
# GeneMark-ES, SNAP, GlimmerHMM and CodingQuarry, aligns the protein evidence
# with Exonerate, and reconciles everything with EvidenceModeler. It picks up
# the transcript evidence from step 05 automatically, from the same directory.
#
# Protein evidence is a curated set of insect proteins plus UniProtKB/Swiss-Prot.
#
# Run from the repository root:
#   conda activate oimp_03_annotation
#   bash 03_annotation/06_funannotate_predict.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
GENOME="${GENOME:-results/03_annotation/repeats/Hifiasm_purged.sorted.fixed.renumbered.fasta.masked}"
OUTDIR="${OUTDIR:-results/03_annotation/funannotate}"
SPECIES="${SPECIES:-Oncideres impluviata}"
PROTEINS="${PROTEINS:-data/insect_proteins_cdhit90.fasta}"   # curated evidence
SWISSPROT="${SWISSPROT:-data/uniprot_sprot.fasta}"
BUSCO_DB="${BUSCO_DB:-insecta}"
PLOIDY="${PLOIDY:-2}"
THREADS="${THREADS:-50}"
# -----------------------------------------------------------------------------

command -v funannotate >/dev/null || { echo "ERROR: funannotate not found" >&2; exit 1; }
[[ -f "$GENOME" ]] || { echo "ERROR: no such file: $GENOME" >&2; exit 1; }

# 1) Predict. --busco_db insecta selects the lineage used to train Augustus.
echo "[1/1] funannotate predict (${THREADS} cpus)"
funannotate predict \
  -i "$GENOME" \
  -s "$SPECIES" \
  --protein_evidence "$PROTEINS" "$SWISSPROT" \
  -o "$OUTDIR" \
  --cpus "$THREADS" \
  --busco_db "$BUSCO_DB" \
  --ploidy "$PLOIDY"

echo "wrote ${OUTDIR}/predict_results/"
