#!/usr/bin/env bash
#
# Orthogroup inference with OrthoFinder, using DIAMOND for the all-versus-all
# searches and MSA mode for orthogroup construction.
#
# The run is timestamped and its results directory recorded in
# orthofinder_out/LATEST_RESULTS_PATH.txt, which every later step reads, so
# re-running does not require any path to be edited downstream.
#
# Run from the repository root:
#   conda activate oimp_04_orthology
#   bash 04_orthology/03_run_orthofinder.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
PROJECT="${PROJECT:-results/04_orthology}"     # section output root
CLEAN="${CLEAN:-${PROJECT}/clean_proteomes}"   # one protein per gene, from step 02
OUTDIR="${OUTDIR:-${PROJECT}/orthofinder_out}" # output directory
THREADS="${THREADS:-64}"                       # CPU threads
ATHREADS="${ATHREADS:-$(( THREADS / 8 > 0 ? THREADS / 8 : 1 ))}"
                                               # threads for the RAM-heavy stage;
                                               # keep well below THREADS
ENV_PREFIX="${ENV_PREFIX:-}"                   # optional conda env prefix to put
                                               # on PATH instead of activating
# -----------------------------------------------------------------------------

# 1) Select the environment. PYTHONPATH and PYTHONHOME are cleared because
#    either one leaks another environment's site-packages into every Python
#    process here, a common cause of import errors. Setting ENV_PREFIX puts an
#    environment's bin/ on PATH directly, which is more robust for batch jobs
#    than `conda activate`.
unset PYTHONPATH PYTHONHOME
if [[ -n "$ENV_PREFIX" ]]; then
  [[ -x "${ENV_PREFIX}/bin/orthofinder" ]] || { echo "Error: no orthofinder at ${ENV_PREFIX}"; exit 1; }
  export PATH="${ENV_PREFIX}/bin:${PATH}"
fi
echo "Using: $(command -v orthofinder)"

# 2) Check the tools. These two are required for the -S diamond -M msa run.
for tool in orthofinder diamond; do
  command -v "$tool" >/dev/null || { echo "Error: $tool not found in PATH"; exit 1; }
done

#    Alternative aligners and tree builders are optional, so warn rather than
#    fail when one is absent.
for tool in famsa mafft fasttree; do
  command -v "$tool" >/dev/null || echo "Note: $tool not on PATH (only needed if you select it)"
done

#    rich is imported by OrthoFinder but missing from the bioconda recipe.
python -c "import rich" 2>/dev/null || {
  echo "Error: python package 'rich' is missing."
  echo "Fix with: ${ENV_PREFIX}/bin/pip install rich"
  exit 1
}

[[ -d "$CLEAN" ]] || { echo "Error: $CLEAN does not exist. Run step 02 first."; exit 1; }
shopt -s nullglob
FA=( "$CLEAN"/*.fa "$CLEAN"/*.faa "$CLEAN"/*.fasta "$CLEAN"/*.pep )
(( ${#FA[@]} >= 2 )) || { echo "Error: found ${#FA[@]} proteomes in $CLEAN (need >=2)"; exit 1; }

RUN="${OUTDIR}/run_$(date +%Y%m%d_%H%M%S)"

echo "== OrthoFinder run =="
echo "  Proteomes:        ${#FA[@]}"
echo "  Search threads:   ${THREADS}  (-t)"
echo "  Analysis threads: ${ATHREADS}  (-a)"
echo "  Output:           ${RUN}"
mkdir -p "$OUTDIR"

# 3) Run. -o must not already exist; OrthoFinder writes straight into it.
orthofinder -f "$CLEAN" \
            -S diamond -M msa \
            -t "$THREADS" -a "$ATHREADS" \
            -o "$RUN"

# 4) Locate the results. With -o they land directly in $RUN, but some versions
#    still nest a Results_* subfolder, so handle both.
RES="$RUN"
if [[ ! -d "${RES}/Orthogroups" ]]; then
  CANDIDATE=$(find "$RUN" -maxdepth 1 -type d -name 'Results_*' | sort | tail -n1)
  [[ -n "$CANDIDATE" ]] && RES="$CANDIDATE"
fi

if [[ ! -d "${RES}/Orthogroups" ]]; then
  echo "Error: could not locate an Orthogroups/ directory under $RUN" >&2
  exit 1
fi

echo "$RES" > "${OUTDIR}/LATEST_RESULTS_PATH.txt"
echo "Latest Results: $RES"

# 5) Summarise, and record the results path for every later step to read.
OG="${RES}/Orthogroups/Orthogroups.tsv"
SCO="${RES}/Orthogroups/Orthogroups_SingleCopyOrthologues.txt"
SCODIR="${RES}/Single_Copy_Orthologue_Sequences"

if [[ -f "$OG" ]]; then
  nSpecies=$(head -n1 "$OG" | awk -F'\t' '{print NF-1}')
  nOG=$(( $(wc -l < "$OG") - 1 ))
  nSCO=$( [[ -f "$SCO" ]] && wc -l < "$SCO" || echo 0 )
  SCOFILES=( "$SCODIR"/*.fa )
  nSCOseq=${#SCOFILES[@]}
  echo "Summary:"
  echo "  Species in matrix:  $nSpecies"
  echo "  Total orthogroups:  $nOG"
  echo "  Single-copy OG IDs: $nSCO"
  echo "  SCO sequence files: $nSCOseq"
else
  echo "Warning: could not find Orthogroups.tsv in $RES"
fi

echo "Done."
