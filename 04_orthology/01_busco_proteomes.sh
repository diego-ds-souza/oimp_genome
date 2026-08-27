#!/usr/bin/env bash
#
# Completeness of every reference proteome, scored with compleasm against the
# insecta_odb10 BUSCO marker set (Table S2).
#
# This runs on the proteomes as downloaded, BEFORE isoform reduction, so the
# duplicated-marker figures reflect the annotations as published rather than
# the collapsed sets used for orthology.
#
# Completed accessions are skipped on a re-run; set FORCE=1 to redo them.
#
# Run from the repository root:
#   conda activate oimp_04_orthology
#   bash 04_orthology/01_busco_proteomes.sh
#
set -uo pipefail

# ------------------------------- settings ------------------------------------
INDIR="${INDIR:-data/proteomes}"                 # NCBI proteomes, as downloaded
OUTDIR="${OUTDIR:-results/04_orthology/compleasm}"
LINEAGE="${LINEAGE:-insecta}"
ODB="${ODB:-odb10}"
THREADS="${THREADS:-16}"
CONDA_ENV="${CONDA_ENV:-}"      # empty: use the already-active environment
LIBRARY_PATH="${LIBRARY_PATH:-}"
FORCE="${FORCE:-0}"                              # 1: re-run completed accessions
OIMP_FAA="${OIMP_FAA:-}"                         # O. impluviata proteins, so the
                                                 # table is complete; empty skips
OIMP_NAME="${OIMP_NAME:-Oncideres_impluviata}"
# -----------------------------------------------------------------------------

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: THREADS must be a positive integer." >&2; exit 1; }
[[ -d "$INDIR" ]] || { echo "ERROR: input directory not found: $INDIR" >&2; exit 1; }

# 1) Make conda available, if an environment was named.
# shellcheck disable=SC1091
if [[ -n "${CONDA_EXE:-}" ]]; then
    source "$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh"
elif [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [[ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
else
    echo "NOTE: conda.sh not found; using the environment already on PATH." >&2
fi

if [[ -n "$CONDA_ENV" ]]; then
    conda activate "$CONDA_ENV" || {
        echo "ERROR: could not activate env '$CONDA_ENV'." >&2
        exit 1
    }
fi
command -v compleasm >/dev/null || {
    echo "ERROR: compleasm not found; activate the section environment first." >&2
    exit 1
}
compleasm protein --help >/dev/null 2>&1 || {
    echo "ERROR: this compleasm installation does not support protein mode." >&2
    exit 1
}

echo "compleasm: $(compleasm --version 2>&1 | head -1)"
echo "lineage : ${LINEAGE}_${ODB}"
echo "threads : $THREADS"
echo "input   : $INDIR"
echo "output  : $OUTDIR"
echo

mkdir -p "$OUTDIR/logs"

COMPLEASM_OPTS=(-l "$LINEAGE" --odb "$ODB" -t "$THREADS")
[[ -n "$LIBRARY_PATH" ]] && COMPLEASM_OPTS+=(-L "$LIBRARY_PATH")

is_complete() {
    [[ -s "$1" ]] && grep -Eq '^N:[[:space:]]*[0-9]+' "$1"
}

# 2) Collect the proteomes to score.
# One entry per proteome: "<name>\t<path to protein.faa>".
# Only files named exactly protein.faa are used, preserving the original
# accession-level selection and ignoring cd-hit/genomic/ncbi_dataset FASTAs.
declare -a NAMES=() PATHS=()
while IFS= read -r faa; do
    acc="$(basename "$(dirname "$faa")")"
    NAMES+=("$acc")
    PATHS+=("$faa")
done < <(find "$INDIR" -mindepth 2 -maxdepth 2 -type f -name 'protein.faa' | sort)

if [[ -n "$OIMP_FAA" ]]; then
    if [[ -f "$OIMP_FAA" ]]; then
        NAMES+=("$OIMP_NAME")
        PATHS+=("$OIMP_FAA")
    else
        echo "WARNING: OIMP_FAA is set but '$OIMP_FAA' does not exist; skipping." >&2
    fi
fi

TOTAL="${#NAMES[@]}"
[[ "$TOTAL" -gt 0 ]] || {
    echo "ERROR: no protein.faa found exactly two levels below $INDIR" >&2
    exit 1
}
echo "found $TOTAL proteomes"
echo

# 3) Run compleasm once per proteome, skipping accessions already done.
declare -a FAILED=()
i=0
for idx in "${!NAMES[@]}"; do
    acc="${NAMES[$idx]}"
    faa="${PATHS[$idx]}"
    run_dir="$OUTDIR/$acc"
    summary="$run_dir/summary.txt"
    run_log="$OUTDIR/logs/$acc.log"
    i=$((i + 1))
    printf '[%2d/%2d] %-20s ' "$i" "$TOTAL" "$acc"

    if is_complete "$summary" && [[ "$FORCE" != "1" ]]; then
        echo "already done, skipping"
        continue
    fi
    [[ "$FORCE" == "1" && -d "$run_dir" ]] && rm -rf -- "${run_dir:?}"

    start=$SECONDS
    if compleasm protein \
        -p "$faa" \
        -o "$run_dir" \
        "${COMPLEASM_OPTS[@]}" \
        >"$run_log" 2>&1 \
        && is_complete "$summary"
    then
        line="$(grep -E '^[SDFM]:' "$summary" | tr '\n' ' ')"
        printf 'done in %4ds   %s\n' "$((SECONDS - start))" "${line:-(summary not parsed)}"
    else
        echo "FAILED  (see logs/$acc.log)"
        FAILED+=("$acc")
    fi
done

# Compleasm downloads a missing lineage unless it is already present in its
# library. For offline use, set LIBRARY_PATH to a populated compleasm library.

# 4) Collate every run into the single table used for Table S2.
TABLE="$OUTDIR/compleasm_summary_all_proteomes.tsv"
{
    printf 'accession\tn_proteins\tcomplete_pct\tcomplete_single_pct\tcomplete_duplicated_pct\t'
    printf 'fragmented_pct\tmissing_pct\tn_markers\tcomplete\tsingle\tduplicated\tfragmented\tmissing\n'
    for idx in "${!NAMES[@]}"; do
        acc="${NAMES[$idx]}"
        faa="${PATHS[$idx]}"
        summary="$OUTDIR/$acc/summary.txt"
        is_complete "$summary" || continue
        np="$(grep -c '^>' "$faa" || true)"
        awk -v acc="$acc" -v np="$np" '
            /^[SDFM]:/ {
                key=substr($1, 1, 1)
                value=$1
                sub(/^[SDFM]:/, "", value)
                sub(/%,$/, "", value)
                pct[key]=value + 0
                count[key]=$2 + 0
            }
            /^N:/ {
                total=$1
                sub(/^N:/, "", total)
            }
            END {
                complete_pct=pct["S"] + pct["D"]
                complete=count["S"] + count["D"]
                printf "%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\t%d\t%d\t%d\t%d\t%d\n", \
                       acc, np, complete_pct, pct["S"], pct["D"], pct["F"], pct["M"], total, \
                       complete, count["S"], count["D"], count["F"], count["M"]
            }
        ' "$summary"
    done
} > "$TABLE"

echo
echo "collated table: $TABLE"
column -t -s $'\t' "$TABLE" 2>/dev/null || cat "$TABLE"

if [[ "${#FAILED[@]}" -gt 0 ]]; then
    echo
    echo "FAILED (${#FAILED[@]}): ${FAILED[*]}"
    exit 1
fi
