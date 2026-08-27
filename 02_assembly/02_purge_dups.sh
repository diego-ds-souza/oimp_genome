#!/usr/bin/env bash
#
# Removal of haplotypic duplications from the primary assembly with purge_dups,
# following the protocol in the hifiasm supplementary material (Cheng et al.
# 2021, Nat Methods 18:170-175).
#
# purge_dups comes with the conda environment. If you would rather build it from
# source, point PD_BIN at the resulting bin directory:
#   git clone https://github.com/dfguan/purge_dups.git
#   cd purge_dups/src && make
#
# Run from the repository root:
#   conda activate oimp_02_assembly
#   bash 02_assembly/02_purge_dups.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
ASM="${ASM:-results/02_assembly/asm-ccsQ20.p_ctg.fasta}"   # primary assembly
READS="${READS:-data/hifi_reads.fastq.gz}"                 # PacBio HiFi reads
OUTDIR="${OUTDIR:-results/02_assembly/purge_dups}"         # output directory
PD_BIN="${PD_BIN:-}"                                       # purge_dups bin dir,
                                                           # empty = use PATH
THREADS="${THREADS:-8}"                                    # CPU threads
# -----------------------------------------------------------------------------

command -v minimap2 >/dev/null || { echo "ERROR: minimap2 not found" >&2; exit 1; }
if [[ -z "$PD_BIN" ]]; then
  PD_BIN=$(dirname "$(command -v pbcstat 2>/dev/null || echo /nonexistent/x)")
fi
[[ -x "$PD_BIN/pbcstat" ]] || { echo "ERROR: purge_dups binaries not found (set PD_BIN)" >&2; exit 1; }
[[ -f "$ASM"   ]] || { echo "ERROR: no such file: $ASM" >&2; exit 1; }
[[ -f "$READS" ]] || { echo "ERROR: no such file: $READS" >&2; exit 1; }

ASM=$(readlink -f "$ASM"); READS=$(readlink -f "$READS")
mkdir -p "$OUTDIR"; cd "$OUTDIR"

# 1) Read depth: map the HiFi reads back to the primary assembly.
echo "[1/6] mapping reads to the assembly"
minimap2 -I6G -xmap-pb "$ASM" "$READS" -t "$THREADS" > read-aln.paf

# 2) Coverage statistics and cutoffs.
#    calcuts is run with its defaults; the cutoffs it chose for this assembly
#    were  5  18  34  35  73  162  (low, mid and high depth thresholds).
#    Inspect PB.stat/PB.cov.wig and rerun with -m/-u only if the automatic
#    cutoffs do not match the coverage histogram.
echo "[2/6] coverage statistics"
"$PD_BIN/pbcstat" read-aln.paf
"$PD_BIN/calcuts" PB.stat > cutoffs 2> calcuts.log

# 3) Self-alignment of the split assembly.
echo "[3/6] splitting the assembly"
"$PD_BIN/split_fa" "$ASM" > split.fa

echo "[4/6] self-aligning the split assembly"
minimap2 -I6G -xasm5 -DP split.fa split.fa -t "$THREADS" > ctg-aln.paf

# 4) Call and remove duplications.
echo "[5/6] identifying duplications"
"$PD_BIN/purge_dups" -2 -T cutoffs -c PB.base.cov ctg-aln.paf > dups.bed

echo "[6/6] extracting the purged assembly"
"$PD_BIN/get_seqs" dups.bed "$ASM"

# get_seqs writes purged.fa (haploid assembly) and hap.fa (removed haplotigs).
mv purged.fa Hifiasm_purged.fa

echo "wrote ${OUTDIR}/Hifiasm_purged.fa  (haplotigs in hap.fa)"
