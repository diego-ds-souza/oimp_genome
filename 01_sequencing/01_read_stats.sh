#!/usr/bin/env bash
#
# Read-level statistics for the PacBio HiFi reads: read count, total bases,
# N50, quality and GC. Divide the total bases by the genome size for coverage.
#
# Run from the repository root:
#   conda activate oimp_01_sequencing
#   bash 01_sequencing/01_read_stats.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
READS="${READS:-data/hifi_reads.fastq.gz}"     # PacBio HiFi reads (gz is fine)
OUTDIR="${OUTDIR:-results/01_sequencing}"      # output directory
OUT="${OUT:-hifi_seqkit_stats.txt}"            # output filename
# -----------------------------------------------------------------------------

command -v seqkit >/dev/null || { echo "ERROR: seqkit not found" >&2; exit 1; }
[[ -f "$READS" ]] || { echo "ERROR: no such file: $READS" >&2; exit 1; }

mkdir -p "$OUTDIR"

# 1) Full statistics: -a adds N50, Q20/Q30 and GC to the default summary.
echo "[1/1] seqkit stats"
seqkit stats -a "$READS" | tee "${OUTDIR}/${OUT}"

echo "wrote ${OUTDIR}/${OUT}"
