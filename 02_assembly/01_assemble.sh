#!/usr/bin/env bash
#
# Primary contig assembly from PacBio HiFi reads with hifiasm, followed by
# conversion of the primary contig graph to FASTA with any2fasta.
#
# Run from the repository root:
#   conda activate oimp_02_assembly
#   bash 02_assembly/01_assemble.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
READS="${READS:-data/hifi_reads.fastq.gz}"     # PacBio HiFi reads (gz is fine)
OUTDIR="${OUTDIR:-results/02_assembly}"        # output directory
PREFIX="${PREFIX:-asm-ccsQ20}"                 # hifiasm output prefix
THREADS="${THREADS:-8}"                        # CPU threads
# -----------------------------------------------------------------------------

for t in hifiasm any2fasta; do
  command -v "$t" >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }
done
[[ -f "$READS" ]] || { echo "ERROR: no such file: $READS" >&2; exit 1; }

mkdir -p "$OUTDIR"
cd "$OUTDIR"
READS_ABS=$(cd "$OLDPWD" && readlink -f "$READS")

# 1) Assemble. --primary keeps the primary and alternate contig graphs; the
#    primary graph (.p_ctg.gfa) is what is carried forward.
echo "[1/2] hifiasm (${THREADS} threads)"
hifiasm -o "$PREFIX" -t "$THREADS" --primary "$READS_ABS"

# 2) Convert the primary contig graph to FASTA.
echo "[2/2] any2fasta"
any2fasta "${PREFIX}.p_ctg.gfa" > "${PREFIX}.p_ctg.fasta"

echo "wrote ${OUTDIR}/${PREFIX}.p_ctg.fasta"
