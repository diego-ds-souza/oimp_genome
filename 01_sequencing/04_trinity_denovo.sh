#!/usr/bin/env bash
#
# Sample-specific de novo transcriptome assemblies with Trinity, from the
# trimmed pairs written by step 03.
#
# These assemblies were supplied to Funannotate together with the genome-guided
# assemblies it builds itself (03_annotation/05_funannotate_train.sh) as the
# transcript evidence behind exon-intron boundaries and UTR refinement.
#
# Trinity is run with its default parameters; only the read pairs, the thread
# count and the memory ceiling are set here.
#
# Run from the repository root:
#   conda activate oimp_01_sequencing
#   bash 01_sequencing/04_trinity_denovo.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
TRIMDIR="${TRIMDIR:-data/rnaseq/trimmed}"      # output of step 03
OUTDIR="${OUTDIR:-results/01_sequencing/trinity}"
THREADS="${THREADS:-16}"                       # CPU threads
MAX_MEMORY="${MAX_MEMORY:-200G}"               # Trinity --max_memory

# The nine libraries of Table S1, as in step 03.
SAMPLES="${SAMPLES:-AF1An AF1Th AF2An AF2Th AM1An AM2An AM2Th L1He L1Mg}"
# -----------------------------------------------------------------------------

command -v Trinity >/dev/null || { echo "ERROR: Trinity not found" >&2; exit 1; }

mkdir -p "$OUTDIR"

# 1) Assemble each library on its own. Trinity requires an output directory
#    whose name contains "trinity"; the assembly is written to Trinity.fasta
#    inside it.
n=0
for s in $SAMPLES; do
  n=$((n + 1))
  echo "[${n}/$(echo $SAMPLES | wc -w)] Trinity ${s}"

  left="${TRIMDIR}/${s}/${s}_R1_fastp.fastq"
  right="${TRIMDIR}/${s}/${s}_R2_fastp.fastq"
  for f in "$left" "$right"; do
    [[ -f "$f" ]] || { echo "ERROR: no such file: $f" >&2; exit 1; }
  done

  Trinity \
    --seqType fq \
    --left "$left" \
    --right "$right" \
    --CPU "$THREADS" \
    --max_memory "$MAX_MEMORY" \
    --output "${OUTDIR}/${s}_trinity"
done

echo "wrote ${OUTDIR}/<sample>_trinity/Trinity.fasta"
