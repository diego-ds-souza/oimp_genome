#!/usr/bin/env bash
#
# Assembly statistics (QUAST) and gene-space completeness (BUSCO), run on both
# the primary assembly and the purged assembly so the effect of haplotig
# removal can be seen. The purged assembly is the one reported in the paper.
#
# Run from the repository root:
#   conda activate oimp_02_assembly
#   bash 02_assembly/03_assembly_qc.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
PRIMARY="${PRIMARY:-results/02_assembly/asm-ccsQ20.p_ctg.fasta}"      # before purging
PURGED="${PURGED:-results/02_assembly/purge_dups/Hifiasm_purged.fa}"  # after purging
OUTDIR="${OUTDIR:-results/02_assembly/qc}"                            # output directory
LINEAGE="${LINEAGE:-insecta_odb10}"                                   # BUSCO lineage
THREADS="${THREADS:-8}"                                               # CPU threads
# -----------------------------------------------------------------------------

for t in quast busco; do
  command -v "$t" >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }
done
for f in "$PRIMARY" "$PURGED"; do
  [[ -f "$f" ]] || { echo "ERROR: no such file: $f" >&2; exit 1; }
done

PRIMARY=$(readlink -f "$PRIMARY"); PURGED=$(readlink -f "$PURGED")
mkdir -p "$OUTDIR"; cd "$OUTDIR"

# 1) Contiguity: length, contig count, N50, GC.
echo "[1/3] QUAST"
quast "$PRIMARY" "$PURGED" -o quast -t "$THREADS"

# 2) Completeness before purging.
echo "[2/3] BUSCO, primary assembly"
busco --in "$PRIMARY" --out busco_primary \
      --lineage_dataset "$LINEAGE" --mode genome --cpu "$THREADS"

# 3) Completeness after purging. Duplicated BUSCOs should drop sharply here;
#    if they do not, the purge_dups cutoffs need revisiting.
echo "[3/3] BUSCO, purged assembly"
busco --in "$PURGED" --out busco_purged \
      --lineage_dataset "$LINEAGE" --mode genome --cpu "$THREADS"

echo "wrote ${OUTDIR}/quast/report.txt"
echo "      ${OUTDIR}/busco_purged/short_summary*.txt"
