#!/usr/bin/env bash
#
# Adapter and quality trimming of the nine RNA-seq libraries with fastp. The
# trimmed pairs are used twice: as transcript evidence for gene prediction
# (section 03) and as the input to quantification (section 07).
#
# Adapters are detected automatically and the quality thresholds are left at
# their defaults. Overlapping pairs are merged, exact duplicates are removed,
# poly-G and poly-X tails are trimmed, and one base is cut from the 3' end of
# read 1. The unmerged pairs are what the later sections read; the merged and
# unpaired reads are written alongside them for reference.
#
# Run from the repository root:
#   conda activate oimp_01_sequencing
#   bash 01_sequencing/03_rnaseq_fastp.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
RAWDIR="${RAWDIR:-data/rnaseq/raw}"            # <sample>/<sample>_R[12].fastq.gz
OUTDIR="${OUTDIR:-data/rnaseq/trimmed}"        # one subdirectory per sample
THREADS="${THREADS:-16}"                       # fastp worker threads

# Sample prefixes: A=adult, L=larva; F/M=female/male; An=antennae, Th=thorax,
# He=head, Mg=midgut. The nine libraries of Table S1.
SAMPLES="${SAMPLES:-AF1An AF1Th AF2An AF2Th AM1An AM2An AM2Th L1He L1Mg}"

TRIM_TAIL1="${TRIM_TAIL1:-1}"                  # bases cut from the 3' end of R1
DUP_ACCURACY="${DUP_ACCURACY:-4}"              # --dup_calc_accuracy, 1 (fast)
                                               # to 6 (most accurate)
# -----------------------------------------------------------------------------

command -v fastp >/dev/null || { echo "ERROR: fastp not found" >&2; exit 1; }

# 1) Check that every library is present before starting, so a missing file
#    does not surface hours into the run.
for s in $SAMPLES; do
  for r in 1 2; do
    f="${RAWDIR}/${s}/${s}_R${r}.fastq.gz"
    [[ -f "$f" ]] || { echo "ERROR: no such file: $f" >&2; exit 1; }
  done
done

# 2) Trim each library in turn. The HTML and JSON reports carry the duplication
#    rate, insert-size distribution and before/after quality profiles.
n=0
for s in $SAMPLES; do
  n=$((n + 1))
  echo "[${n}/$(echo $SAMPLES | wc -w)] fastp ${s}"

  r1_file="${RAWDIR}/${s}/${s}_R1.fastq.gz"
  r2_file="${RAWDIR}/${s}/${s}_R2.fastq.gz"

  sample_dir="${OUTDIR}/${s}"
  mkdir -p "$sample_dir"

  out_r1="${sample_dir}/${s}_R1_fastp.fastq"
  out_r2="${sample_dir}/${s}_R2_fastp.fastq"
  unpaired="${sample_dir}/${s}_unpaired_fastp.fastq"
  merged="${sample_dir}/${s}_merged_fastp.fastq"
  html_report="${sample_dir}/${s}_fastp.html"
  json_report="${sample_dir}/${s}_fastp.json"
  log_file="${sample_dir}/${s}_fastp.log"

  fastp \
    -i "$r1_file" \
    -I "$r2_file" \
    -o "$out_r1" \
    -O "$out_r2" \
    --unpaired1 "$unpaired" \
    --merge \
    --merged_out "$merged" \
    --trim_tail1 "$TRIM_TAIL1" \
    --dedup \
    --dup_calc_accuracy "$DUP_ACCURACY" \
    --trim_poly_g \
    --trim_poly_x \
    --html "$html_report" \
    --json "$json_report" \
    -w "$THREADS" \
    2>&1 | tee -a "$log_file"
done

echo "wrote ${OUTDIR}/<sample>/<sample>_R[12]_fastp.fastq"
