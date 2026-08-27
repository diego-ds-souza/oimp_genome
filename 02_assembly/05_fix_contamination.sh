#!/usr/bin/env bash
#
# Apply the NCBI contamination screen results to the sorted assembly and
# renumber the scaffolds. Two edits were required:
#
#   - scaffold_29 was a contaminant and was removed entirely;
#   - an adaptor match inside scaffold_3 was internal, so that scaffold was
#     split at the adaptor into scaffold_3a and scaffold_3b.
#
# Scaffolds are then renumbered scaffold_1 .. scaffold_N by length, giving the
# names used in the paper and in the GenBank assembly (GCA_060235215.1).
#
# The annotation GFF3 is updated alongside the FASTA if GFF_IN is set, so this
# script is normally run twice: once on the assembly alone, and again after
# annotation (section 03) to bring the GFF into the same coordinate system.
# Features spanning the adaptor cut are dropped and logged.
#
# Run from the repository root:
#   conda activate oimp_02_assembly
#   bash 02_assembly/05_fix_contamination.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
FASTA_IN="${FASTA_IN:-results/02_assembly/Hifiasm_purged.sorted.fasta}"
GFF_IN="${GFF_IN:-}"                  # optional annotation GFF3; empty to skip
OUTDIR="${OUTDIR:-results/02_assembly/fixed}"

CUT_SCAFF="${CUT_SCAFF:-scaffold_3}"  # scaffold carrying the internal adaptor
CUT_START="${CUT_START:-39412980}"    # adaptor interval, 1-based inclusive
CUT_END="${CUT_END:-39413053}"
DROP_SCAFF="${DROP_SCAFF:-scaffold_29}"   # contaminant scaffold, removed
# -----------------------------------------------------------------------------

for t in samtools seqkit awk sed; do
  command -v "$t" >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }
done
[[ -f "$FASTA_IN" ]] || { echo "ERROR: no such file: $FASTA_IN" >&2; exit 1; }
if [[ -n "$GFF_IN" ]]; then
  [[ -f "$GFF_IN" ]] || { echo "ERROR: no such file: $GFF_IN" >&2; exit 1; }
  for t in gffread agat_sq_rename_seqid.pl; do
    command -v "$t" >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }
  done
fi

mkdir -p "$OUTDIR"/{logs,tmp}
LOG="$OUTDIR/logs/run.log"; : > "$LOG"
echo "FASTA: $FASTA_IN" | tee -a "$LOG"
[[ -n "$GFF_IN" ]] && echo "GFF  : $GFF_IN" | tee -a "$LOG"

# ---- 1. baseline lengths ----------------------------------------------------
samtools faidx "$FASTA_IN"
cut -f1,2 "${FASTA_IN}.fai" > "$OUTDIR/logs/baseline_lengths.tsv"

SCAFF_LEN=$(awk -v s="$CUT_SCAFF" '$1==s{print $2}' "${FASTA_IN}.fai")
[[ -n "${SCAFF_LEN:-}" ]] || { echo "ERROR: $CUT_SCAFF not in the FASTA index" >&2; exit 1; }
LEFT_END=$((CUT_START - 1))
RIGHT_START=$((CUT_END + 1))
(( LEFT_END >= 1 && RIGHT_START <= SCAFF_LEN )) || { echo "ERROR: cut out of range" >&2; exit 1; }
echo "cutting $CUT_SCAFF ($SCAFF_LEN bp) at ${CUT_START}-${CUT_END}" | tee -a "$LOG"
echo "  -> left 1..${LEFT_END}, right ${RIGHT_START}..${SCAFF_LEN}" | tee -a "$LOG"

# ---- 2. fixed FASTA: drop the contaminant, split the adaptor scaffold -------
seqkit grep -r -v -p "$DROP_SCAFF" "$FASTA_IN" \
  | seqkit grep -r -v -p "$CUT_SCAFF" > "$OUTDIR/tmp/other.fa"
samtools faidx "$FASTA_IN" "${CUT_SCAFF}:1-${LEFT_END}"              > "$OUTDIR/tmp/left.fa"
samtools faidx "$FASTA_IN" "${CUT_SCAFF}:${RIGHT_START}-${SCAFF_LEN}" > "$OUTDIR/tmp/right.fa"
sed -i 's/\r$//' "$OUTDIR/tmp/"*.fa
seqkit replace -p '.*' -r "${CUT_SCAFF}a" "$OUTDIR/tmp/left.fa"  > "$OUTDIR/tmp/left_ren.fa"
seqkit replace -p '.*' -r "${CUT_SCAFF}b" "$OUTDIR/tmp/right.fa" > "$OUTDIR/tmp/right_ren.fa"

FASTA_FIXED="$OUTDIR/tmp/fixed.fasta"
cat "$OUTDIR/tmp/other.fa" "$OUTDIR/tmp/left_ren.fa" "$OUTDIR/tmp/right_ren.fa" > "$FASTA_FIXED"
printf "\n" >> "$FASTA_FIXED"
samtools faidx "$FASTA_FIXED"

# ---- 3. fixed GFF: split the adaptor scaffold, drop spanning features -------
GFF_FIXED=""
if [[ -n "$GFF_IN" ]]; then
  T="$OUTDIR/tmp"
  awk -v SCAF="$CUT_SCAFF" -v L="$LEFT_END" -v R="$RIGHT_START" -v T="$T" -F'\t' '
    BEGIN{OFS="\t"}
    /^#/ { print > (T"/HDR.gff3"); next }
    {
      if ($1 != SCAF) { print > (T"/OTHER.gff3"); next }
      if ($5 <= L)    { print > (T"/LEFT.gff3");  next }
      if ($4 >= R)    { print > (T"/RIGHT.gff3"); next }
      print > (T"/OVERLAP.gff3")
    }' "$GFF_IN"
  for f in HDR OTHER LEFT RIGHT OVERLAP; do [[ -f "$T/$f.gff3" ]] || : > "$T/$f.gff3"; done

  awk -v S="$CUT_SCAFF" -F'\t' 'BEGIN{OFS="\t"} !/^#/ { $1=S"a" } {print}' \
    "$T/LEFT.gff3" > "$T/LEFT.renamed.gff3"
  awk -v S="$CUT_SCAFF" -v SHIFT="$CUT_END" -F'\t' 'BEGIN{OFS="\t"}
    !/^#/ { $1=S"b"; $4=$4-SHIFT; $5=$5-SHIFT; if($4<1) $4=1 } {print}' \
    "$T/RIGHT.gff3" > "$T/RIGHT.rebased.gff3"

  awk -F'\t' 'BEGIN{OFS="\t"} !/^#/{print $1,$3,$4,$5,$9}' "$T/OVERLAP.gff3" \
    > "$OUTDIR/logs/features_overlapping_cut.tsv"
  N_OV=$(wc -l < "$OUTDIR/logs/features_overlapping_cut.tsv")
  echo "features spanning the cut (dropped): $N_OV" | tee -a "$LOG"

  GFF_FIXED="$T/fixed.gff3"
  cat "$T/HDR.gff3" "$T/OTHER.gff3" "$T/LEFT.renamed.gff3" "$T/RIGHT.rebased.gff3" > "$GFF_FIXED"
  gffread -E "$GFF_FIXED" -g "$FASTA_FIXED" -o /dev/null 2>>"$LOG" && echo "gffread OK" | tee -a "$LOG"
fi

# ---- 4. renumber scaffold_1 .. scaffold_N by length -------------------------
MAP="$OUTDIR/logs/renaming_map.tsv"
cut -f1,2 "${FASTA_FIXED}.fai" | sort -k2,2nr \
  | awk 'BEGIN{OFS="\t"} {print $1,"scaffold_"NR,$2}' > "$MAP"

FASTA_OUT="$OUTDIR/Hifiasm_purged.sorted.fixed.renumbered.fasta"
awk 'NR==FNR{m[$1]=$2; next}
     /^>/{h=substr($0,2); sub(/[ \t].*$/,"",h); print (h in m ? ">"m[h] : $0); next}
     {print}' "$MAP" "$FASTA_FIXED" > "$FASTA_OUT"
samtools faidx "$FASTA_OUT"
cut -f1,2 "${FASTA_OUT}.fai" > "$OUTDIR/logs/final_lengths.tsv"

if [[ -n "$GFF_FIXED" ]]; then
  GFF_OUT="$OUTDIR/$(basename "${GFF_IN%.gff3}").fixed.renumbered.gff3"
  agat_sq_rename_seqid.pl --gff "$GFF_FIXED" --tsv "$MAP" -o "$GFF_OUT" >/dev/null
  gffread -E "$GFF_OUT" -g "$FASTA_OUT" -o /dev/null 2>>"$LOG" && echo "final gffread OK" | tee -a "$LOG"
fi

# ---- 5. changelog -----------------------------------------------------------
{
  echo "Dropped   : $DROP_SCAFF"
  echo "Split     : $CUT_SCAFF at ${CUT_START}-${CUT_END}"
  echo
  echo "Renaming map (old -> new, by length)"
  printf 'old_id\tnew_id\tlength_bp\n'
  cat "$MAP"
} > "$OUTDIR/logs/CHANGELOG.txt"

rm -rf "$OUTDIR/tmp"
echo "wrote $FASTA_OUT"
[[ -n "$GFF_FIXED" ]] && echo "      $GFF_OUT"
echo "      $OUTDIR/logs/CHANGELOG.txt"
