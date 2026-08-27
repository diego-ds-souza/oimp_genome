#!/usr/bin/env bash
#
# Transcript evidence for gene prediction. `funannotate train` assembles the
# RNA-seq de novo with Trinity, aligns it to the genome with GMAP and BLAT,
# runs a genome-guided assembly through HISAT2 and StringTie, and refines the
# alignments with PASA. Its output directory is the input to `predict`.
#
# The nine libraries are the adult antennae and thorax and the larval head and
# midgut samples listed in Table S1.
#
# Run from the repository root:
#   conda activate oimp_03_annotation
#   bash 03_annotation/05_funannotate_train.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
GENOME="${GENOME:-results/03_annotation/repeats/Hifiasm_purged.sorted.fixed.renumbered.fasta.masked}"
RNADIR="${RNADIR:-data/rnaseq/trimmed}"   # trimmed reads from 01_sequencing/03
OUTDIR="${OUTDIR:-results/03_annotation/funannotate}"
SPECIES="${SPECIES:-Oncideres impluviata}"
THREADS="${THREADS:-50}"

# Sample prefixes: A=adult, L=larva; F/M=female/male; An=antennae, Th=thorax,
# He=head, Mg=midgut.
SAMPLES="${SAMPLES:-AF1An AF1Th AF2An AF2Th AM1An AM2An AM2Th L1He L1Mg}"
# -----------------------------------------------------------------------------

command -v funannotate >/dev/null || { echo "ERROR: funannotate not found" >&2; exit 1; }
[[ -f "$GENOME" ]] || { echo "ERROR: no such file: $GENOME" >&2; exit 1; }

# 1) Collect the read pairs in a fixed order.
left=(); right=()
for s in $SAMPLES; do
  l="${RNADIR}/${s}/${s}_R1_fastp.fastq"; r="${RNADIR}/${s}/${s}_R2_fastp.fastq"
  [[ -f "$l" && -f "$r" ]] || { echo "ERROR: missing reads for $s" >&2; exit 1; }
  left+=("$l"); right+=("$r")
done

# 2) Train. --stranded no because the libraries are unstranded.
echo "[1/1] funannotate train (${THREADS} cpus, ${#left[@]} libraries)"
funannotate train \
  -i "$GENOME" \
  -o "$OUTDIR" \
  --left "${left[@]}" \
  --right "${right[@]}" \
  --stranded no \
  --species "$SPECIES" \
  --cpus "$THREADS"

echo "wrote ${OUTDIR}/training/"
