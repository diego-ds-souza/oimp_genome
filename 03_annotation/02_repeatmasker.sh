#!/usr/bin/env bash
#
# Annotate and soft-mask repeats with RepeatMasker using the custom library
# from step 01, then summarise the Kimura divergence of every repeat family.
#
# -a writes the alignment file that calcDivergenceFromAlign.pl needs; without
# it there is no divergence summary and no repeat landscape.
#
# Run from the repository root:
#   conda activate oimp_03_annotation
#   bash 03_annotation/02_repeatmasker.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
GENOME="${GENOME:-results/02_assembly/fixed/Hifiasm_purged.sorted.fixed.renumbered.fasta}"
LIB="${LIB:-results/03_annotation/repeats/custom_db.fa}"
OUTDIR="${OUTDIR:-results/03_annotation/repeats}"
NAME="${NAME:-Oimp}"
PA_JOBS="${PA_JOBS:-10}"                 # parallel jobs; rmblast uses
                                         # 4 CPUs per job, so 10 = 40 cores
# -----------------------------------------------------------------------------

command -v RepeatMasker >/dev/null || { echo "ERROR: RepeatMasker not found" >&2; exit 1; }
for f in "$GENOME" "$LIB"; do
  [[ -f "$f" ]] || { echo "ERROR: no such file: $f" >&2; exit 1; }
done

GENOME=$(readlink -f "$GENOME"); LIB=$(readlink -f "$LIB")
mkdir -p "$OUTDIR"
base=$(basename "$GENOME")

# 1) Mask. -s is the slow, most sensitive search; -xsmall soft-masks
#    (lowercase) as funannotate expects; -gff writes the repeat annotation;
#    -a keeps the alignments that step 2 needs.
echo "[1/2] RepeatMasker (${PA_JOBS} parallel jobs)"
RepeatMasker -e rmblast -pa "$PA_JOBS" -s -a -xsmall -gff \
             -dir "$OUTDIR" -lib "$LIB" "$GENOME"

# 2) Kimura divergence per family, corrected for CpG sites. The .divsum holds
#    both the per-family table and the class-by-divergence matrix used by the
#    repeat landscape.
echo "[2/2] calcDivergenceFromAlign.pl"
calcDivergenceFromAlign.pl -s "${OUTDIR}/${NAME}.divsum" "${OUTDIR}/${base}.align"

echo "wrote ${OUTDIR}/${base}.masked"
echo "      ${OUTDIR}/${base}.tbl"
echo "      ${OUTDIR}/${NAME}.divsum"
