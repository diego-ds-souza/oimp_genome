#!/usr/bin/env bash
#
# Quantify every RNA-seq library against the Salmon index from step 02.
#
# Quantification uses automatic library-type detection and Salmon's
# sequence-specific, GC and positional bias corrections. Samples are read from
# samples.tsv, which names the fastq pairs and the stage, sex and tissue of
# each library (Table S1).
#
# Run from the repository root:
#   conda activate oimp_07_transcriptomics
#   bash 07_transcriptomics/03_salmon_quant.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
INDEX="${INDEX:-results/07_transcriptomics/salmon_index}"   # from step 02
SAMPLES="${SAMPLES:-07_transcriptomics/samples.tsv}"        # library metadata
OUTDIR="${OUTDIR:-results/07_transcriptomics/quant}"        # output directory
THREADS="${THREADS:-8}"                                     # CPU threads
# -----------------------------------------------------------------------------

LOGDIR="$(dirname "$OUTDIR")/logs"
mkdir -p "$OUTDIR" "$LOGDIR"

STAMP=$(date +"%Y%m%d_%H%M%S")
LOG="${LOGDIR}/03_salmon_quant_${STAMP}.log"

echo "=============================================================================" | tee "$LOG"
echo "Salmon Quantification" | tee -a "$LOG"
echo "=============================================================================" | tee -a "$LOG"
echo "Index:   $INDEX" | tee -a "$LOG"
echo "Samples: $SAMPLES" | tee -a "$LOG"
echo "Output:  $OUTDIR" | tee -a "$LOG"
echo "Started: $(date)" | tee -a "$LOG"
echo "" | tee -a "$LOG"

# 1) Verify inputs.
if [[ ! -d "$INDEX" ]]; then
    echo "[ERROR] Salmon index not found: $INDEX" | tee -a "$LOG"
    echo "        Run 02_salmon_index.sh first." | tee -a "$LOG"
    exit 1
fi

if [[ ! -f "$SAMPLES" ]]; then
    echo "[ERROR] Sample metadata not found: $SAMPLES" | tee -a "$LOG"
    exit 1
fi

# Check Salmon is installed
if ! command -v salmon &> /dev/null; then
    echo "[ERROR] Salmon not found. Install with: conda install -c bioconda salmon" | tee -a "$LOG"
    exit 1
fi

# Determine number of threads
if command -v nproc &> /dev/null; then
    THREADS=$(nproc)
elif command -v sysctl &> /dev/null; then
    THREADS=$(sysctl -n hw.logicalcpu)
else
    THREADS=4
fi

echo "[INFO] Using $THREADS threads per sample" | tee -a "$LOG"
echo "" | tee -a "$LOG"

# 2) Process each sample.
n_samples=0
n_success=0
n_skipped=0
n_failed=0

# Read samples.tsv (skip header)
# Format: sample_id  stage  sex  structure  group  batch  fastq1  fastq2
tail -n +2 "$SAMPLES" | while IFS=$'\t' read -r sid stage sex structure group batch r1 r2; do
    
    n_samples=$((n_samples + 1))
    
    echo "----------------------------------------" | tee -a "$LOG"
    echo "[Sample $n_samples] $sid" | tee -a "$LOG"
    echo "  Stage:     $stage" | tee -a "$LOG"
    echo "  Sex:       $sex" | tee -a "$LOG"
    echo "  Structure: $structure" | tee -a "$LOG"
    echo "  Group:     $group" | tee -a "$LOG"
    
    sample_out="${OUTDIR}/${sid}"
    
    # Skip if already quantified
    if [[ -s "${sample_out}/quant.sf" ]]; then
        echo "  [SKIP] quant.sf already exists" | tee -a "$LOG"
        n_skipped=$((n_skipped + 1))
        continue
    fi
    
    # Resolve FASTQ paths. Relative paths are taken from the repository root,
    # which is where every script in this repository is run from.
    fq1="$r1"
    fq2="$r2"
    
    # Check FASTQ files exist
    if [[ ! -f "$fq1" ]]; then
        echo "  [ERROR] R1 not found: $fq1" | tee -a "$LOG"
        n_failed=$((n_failed + 1))
        continue
    fi
    if [[ ! -f "$fq2" ]]; then
        echo "  [ERROR] R2 not found: $fq2" | tee -a "$LOG"
        n_failed=$((n_failed + 1))
        continue
    fi
    
    echo "  R1: $fq1" | tee -a "$LOG"
    echo "  R2: $fq2" | tee -a "$LOG"
    echo "  [RUNNING] Salmon quant..." | tee -a "$LOG"
    
    # Run Salmon quantification
    if salmon quant \
        -i "$INDEX" \
        -l A \
        -1 "$fq1" \
        -2 "$fq2" \
        -p "$THREADS" \
        --validateMappings \
        --seqBias \
        --gcBias \
        --posBias \
        -o "$sample_out" \
        2>> "$LOG"; then
        
        echo "  [SUCCESS] Quantification complete" | tee -a "$LOG"
        n_success=$((n_success + 1))
        
        # Report mapping rate
        if [[ -f "${sample_out}/aux_info/meta_info.json" ]]; then
            map_rate=$(grep -o '"percent_mapped": [0-9.]*' "${sample_out}/aux_info/meta_info.json" | cut -d' ' -f2)
            echo "  Mapping rate: ${map_rate}%" | tee -a "$LOG"
        fi
    else
        echo "  [FAILED] Salmon quant failed" | tee -a "$LOG"
        n_failed=$((n_failed + 1))
    fi
    
done

# 3) Summary.
echo "" | tee -a "$LOG"
echo "=============================================================================" | tee -a "$LOG"
echo "SUMMARY" | tee -a "$LOG"
echo "=============================================================================" | tee -a "$LOG"
echo "Samples processed: $n_samples" | tee -a "$LOG"
echo "  Successful: $n_success" | tee -a "$LOG"
echo "  Skipped:    $n_skipped" | tee -a "$LOG"
echo "  Failed:     $n_failed" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "Output directory: $OUTDIR" | tee -a "$LOG"
echo "Completed: $(date)" | tee -a "$LOG"

# List output
echo "" | tee -a "$LOG"
echo "Quantification files:" | tee -a "$LOG"
ls -la "$OUTDIR"/*/quant.sf 2>/dev/null | tee -a "$LOG" || echo "  No quant.sf files found" | tee -a "$LOG"
