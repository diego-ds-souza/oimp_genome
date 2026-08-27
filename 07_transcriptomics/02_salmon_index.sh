#!/usr/bin/env bash
#
# Build the Salmon index used for quantification, with the genome supplied as
# decoy sequence. Decoy-aware indexing stops reads from genomic regions that
# resemble a transcript being assigned to that transcript, which otherwise
# inflates the counts of the affected genes.
#
# Needs roughly 16 GB of memory.
#
# Run from the repository root:
#   conda activate oimp_07_transcriptomics
#   bash 07_transcriptomics/02_salmon_index.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
GENOME="${GENOME:-results/02_assembly/fixed/Hifiasm_purged.sorted.fixed.renumbered.fasta}"
TRANSCRIPTS="${TRANSCRIPTS:-results/03_annotation/curated/Oncideres_impluviata.mrna.fa}"
OUTDIR="${OUTDIR:-results/07_transcriptomics}"   # output directory
INDEX="${INDEX:-${OUTDIR}/salmon_index}"         # index written here
THREADS="${THREADS:-8}"                          # CPU threads
# -----------------------------------------------------------------------------

REFDIR="$OUTDIR"
LOGDIR="${OUTDIR}/logs"
mkdir -p "$REFDIR" "$LOGDIR"

STAMP=$(date +"%Y%m%d_%H%M%S")
LOG="${LOGDIR}/02_salmon_index_${STAMP}.log"

echo "=============================================================================" | tee "$LOG"
echo "Build Salmon Index with Decoys" | tee -a "$LOG"
echo "=============================================================================" | tee -a "$LOG"
echo "Transcripts: $TRANSCRIPTS" | tee -a "$LOG"
echo "Genome:      $GENOME" | tee -a "$LOG"
echo "Index:       $INDEX" | tee -a "$LOG"
echo "Started:     $(date)" | tee -a "$LOG"
echo "" | tee -a "$LOG"

# 1) Verify inputs.
if [[ ! -f "$TRANSCRIPTS" ]]; then
    echo "[ERROR] Transcript FASTA not found: $TRANSCRIPTS" | tee -a "$LOG"
    exit 1
fi

if [[ ! -f "$GENOME" ]]; then
    echo "[ERROR] Genome FASTA not found: $GENOME" | tee -a "$LOG"
    exit 1
fi

# Check Salmon is installed
if ! command -v salmon &> /dev/null; then
    echo "[ERROR] Salmon not found. Install with: conda install -c bioconda salmon" | tee -a "$LOG"
    exit 1
fi

echo "[INFO] Salmon version:" | tee -a "$LOG"
salmon --version 2>&1 | tee -a "$LOG"
echo "" | tee -a "$LOG"

# 2) Create decoy list (scaffold/contig names from genome).
echo "[Step 1] Creating decoy list from genome..." | tee -a "$LOG"

DECOYS="${REFDIR}/decoys.txt"
grep "^>" "$GENOME" | sed 's/^>//' | cut -d ' ' -f1 > "$DECOYS"

n_decoys=$(wc -l < "$DECOYS" | tr -d ' ')
echo "  Decoy sequences: $n_decoys" | tee -a "$LOG"

# 3) Create gentrome (transcripts + genome).
echo "[Step 2] Creating gentrome (transcripts + genome)..." | tee -a "$LOG"

GENTROME="${REFDIR}/gentrome.fa"

# Concatenate transcripts first, then genome
cat "$TRANSCRIPTS" "$GENOME" > "$GENTROME"

# Get sizes
tx_size=$(du -h "$TRANSCRIPTS" | cut -f1)
genome_size=$(du -h "$GENOME" | cut -f1)
gentrome_size=$(du -h "$GENTROME" | cut -f1)

echo "  Transcripts size: $tx_size" | tee -a "$LOG"
echo "  Genome size:      $genome_size" | tee -a "$LOG"
echo "  Gentrome size:    $gentrome_size" | tee -a "$LOG"

# 4) Build salmon index.
echo "[Step 3] Building Salmon index (this may take 10-30 minutes)..." | tee -a "$LOG"

# Determine number of threads
if command -v nproc &> /dev/null; then
    THREADS=$(nproc)
elif command -v sysctl &> /dev/null; then
    THREADS=$(sysctl -n hw.logicalcpu)
else
    THREADS=4
fi

echo "  Using $THREADS threads" | tee -a "$LOG"

salmon index \
    -t "$GENTROME" \
    -d "$DECOYS" \
    -i "$INDEX" \
    -k 31 \
    --threads "$THREADS" \
    2>&1 | tee -a "$LOG"

# 5) Verify output.
echo "" | tee -a "$LOG"
if [[ -d "$INDEX" && -f "$INDEX/info.json" ]]; then
    echo "[SUCCESS] Salmon index created successfully!" | tee -a "$LOG"
    echo "" | tee -a "$LOG"
    echo "Index contents:" | tee -a "$LOG"
    ls -lh "$INDEX" | tee -a "$LOG"
    echo "" | tee -a "$LOG"
    echo "Index info:" | tee -a "$LOG"
    cat "$INDEX/info.json" | tee -a "$LOG"
else
    echo "[ERROR] Salmon index creation failed!" | tee -a "$LOG"
    exit 1
fi

# 6) Cleanup (optional - comment out to keep intermediate files).
# echo "[Cleanup] Removing intermediate gentrome file..." | tee -a "$LOG"
# rm -f "$GENTROME"

echo "" | tee -a "$LOG"
echo "=============================================================================" | tee -a "$LOG"
echo "DONE!" | tee -a "$LOG"
echo "=============================================================================" | tee -a "$LOG"
echo "Index location: $INDEX" | tee -a "$LOG"
echo "Completed: $(date)" | tee -a "$LOG"
