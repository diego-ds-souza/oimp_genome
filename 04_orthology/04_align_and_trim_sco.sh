#!/usr/bin/env bash
#
# Align each single-copy orthologue with MAFFT and trim the alignment with
# trimAl, giving the per-gene alignments that step 05 concatenates.
#
# Run from the repository root:
#   conda activate oimp_04_orthology
#   bash 04_orthology/04_align_and_trim_sco.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
PROJECT="${PROJECT:-results/04_orthology}"        # section output root
OF_BASE="${OF_BASE:-${PROJECT}/orthofinder_out}"  # OrthoFinder output, from step 03
ALIGN_DIR="${ALIGN_DIR:-${PROJECT}/sco_align}"    # MAFFT alignments
TRIM_DIR="${TRIM_DIR:-${PROJECT}/trimmed}"        # trimAl output
THREADS="${THREADS:-16}"                          # CPU threads
# -----------------------------------------------------------------------------

command -v mafft  >/dev/null || { echo "Error: mafft not found in PATH"; exit 1; }
command -v trimal >/dev/null || { echo "Error: trimal not found in PATH"; exit 1; }

# 1) Locate the OrthoFinder results: the path step 03 recorded if it exists,
#    otherwise the newest Results_* directory found.
if [[ -f "${OF_BASE}/LATEST_RESULTS_PATH.txt" ]]; then
  RES="$(cat "${OF_BASE}/LATEST_RESULTS_PATH.txt")"
else
  RES=$(ls -d ${OF_BASE}/run_*/Results_* 2>/dev/null | sort | tail -n1 || true)
  if [[ -z "${RES}" ]]; then
    RES=$(ls -d ${OF_BASE}/Results_* 2>/dev/null | sort | tail -n1 || true)
  fi
fi

[[ -n "${RES:-}" && -d "$RES" ]] || { echo "Error: could not locate a Results_* folder under ${OF_BASE}"; exit 1; }

SCO_SRC="${RES}/Single_Copy_Orthologue_Sequences"
[[ -d "$SCO_SRC" ]] || { echo "Error: ${SCO_SRC} not found. Did OrthoFinder finish?"; exit 1; }

echo "Using OrthoFinder results: $RES"
echo "SCO source:               $SCO_SRC"

mkdir -p "$ALIGN_DIR" "$TRIM_DIR"

# 2) Copy the single-copy orthologue sequences across. -n leaves any existing
#    file alone, so a re-run does not clobber finished alignments.
shopt -s nullglob
cp -n "${SCO_SRC}"/*.fa "$ALIGN_DIR" || true

# 3) Align. Orthogroups that already have an .aln.fa are skipped, so an
#    interrupted run resumes where it stopped.
echo "Aligning with MAFFT using ${THREADS} threads..."
n_aligned=0
for f in "${ALIGN_DIR}"/*.fa; do
  [[ -f "${f%.fa}.aln.fa" ]] && continue
  mafft --thread "${THREADS}" --auto "$f" > "${f%.fa}.aln.fa"
  ((n_aligned++)) || true
done
echo "MAFFT: created ${n_aligned} new alignments (others already present)."

# 4) Trim with the automated1 heuristic, skipping alignments already done.
echo "Trimming with trimAl..."
n_trimmed=0
for f in "${ALIGN_DIR}"/*.aln.fa; do
  out="${TRIM_DIR}/$(basename "$f")"
  [[ -f "$out" ]] && continue
  trimal -in "$f" -out "$out" -automated1
  ((n_trimmed++)) || true
done
echo "trimAl: created ${n_trimmed} new trimmed alignments."

# 5) Report the counts at each stage; these three should agree.
n_sco=$(ls -1 "${SCO_SRC}"/*.fa 2>/dev/null | wc -l | tr -d ' ')
n_aln=$(ls -1 "${ALIGN_DIR}"/*.aln.fa 2>/dev/null | wc -l | tr -d ' ')
n_tri=$(ls -1 "${TRIM_DIR}"/*.aln.fa 2>/dev/null | wc -l | tr -d ' ')
echo "Summary:"
echo "  SCO files:       ${n_sco}"
echo "  MAFFT alignments:${n_aln}"
echo "  Trimmed files:   ${n_tri}"

echo "Done."
