#!/usr/bin/env bash
#
# k-mer counting for genome-size estimation. The histogram written here is the
# input to GenomeScope (see README.md).
#
# Run from the repository root:
#   conda activate oimp_01_sequencing
#   bash 01_sequencing/02_kmer_count.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
READS="${READS:-data/hifi_reads.fastq}"        # PacBio HiFi reads, uncompressed
OUTDIR="${OUTDIR:-results/01_sequencing}"      # output directory
K="${K:-21}"                                   # k-mer length (21 for GenomeScope)
THREADS="${THREADS:-10}"                       # CPU threads
HASH="${HASH:-1000000000}"                     # hash size; raise this if
                                               # jellyfish reports it is full
# -----------------------------------------------------------------------------

command -v jellyfish >/dev/null || { echo "ERROR: jellyfish not found" >&2; exit 1; }
[[ -f "$READS" ]] || { echo "ERROR: no such file: $READS" >&2; exit 1; }

mkdir -p "$OUTDIR"
name=$(basename "$READS"); name=${name%%.*}
prefix="${OUTDIR}/${name}_k${K}"

# 1) Count k-mers. -C counts canonical k-mers, which GenomeScope requires.
echo "[1/2] counting ${K}-mers (canonical, ${THREADS} threads, hash ${HASH})"
jellyfish count -C -m "$K" -s "$HASH" -t "$THREADS" -o "${prefix}.jf" "$READS"

# 2) Collapse the counts into the coverage histogram GenomeScope reads.
echo "[2/2] building histogram"
jellyfish histo -t "$THREADS" "${prefix}.jf" > "${prefix}.histo"

echo "wrote ${prefix}.histo  -> upload to GenomeScope"
