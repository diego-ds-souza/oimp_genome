#!/usr/bin/env bash
#
# Reduce every downloaded proteome to one protein per gene, keeping the longest
# isoform. NCBI protein FASTA headers carry no gene identifier, so each protein
# is mapped back to its parent gene through the genomic GFF3; see
# isoform_collapse/README.md for why a header-based shortcut fails silently.
#
# Comparative methods downstream - OrthoFinder, CAFE, gene family counting -
# all assume one sequence per gene, so this step decides their input.
#
# Run from the repository root:
#   conda activate oimp_04_orthology
#   bash 04_orthology/02_collapse_isoforms.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
INDIR="${INDIR:-data/proteomes}"                        # from step 00
OIMP="${OIMP:-results/03_annotation/curated/Oncideres_impluviata.proteins.clean.fa}"
OUTDIR="${OUTDIR:-results/04_orthology/clean_proteomes}"
MANIFEST="${MANIFEST:-results/04_orthology/manifest.tsv}"
MIN_LENGTH="${MIN_LENGTH:-30}"                          # drop shorter proteins
# -----------------------------------------------------------------------------

HERE="04_orthology/isoform_collapse"
[[ -d "$INDIR" ]] || { echo "ERROR: no such directory: $INDIR" >&2; exit 1; }
mkdir -p "$(dirname "$MANIFEST")" "$OUTDIR"

# 1) Build the manifest: one line per genome, giving the protein FASTA and the
#    GFF3 to read gene identifiers from.
echo "[1/3] building the manifest"
python3 "${HERE}/make_manifest.py" "$INDIR" > "$MANIFEST"

# 2) O. impluviata is already one protein per gene, so it is listed as
#    pre-collapsed ("-" in the GFF3 field) and copied through unchanged.
if [[ -f "$OIMP" ]]; then
  printf 'O_impluviata\tO_impluviata\t%s\t-\n' "$OIMP" >> "$MANIFEST"
fi

# 3) Collapse, then verify. The check in step 3 is the one that matters: it
#    catches a reduction that silently did nothing.
echo "[2/3] collapsing isoforms"
bash "${HERE}/run_collapse_batch.sh" "$MANIFEST" "$OUTDIR" --min-length "$MIN_LENGTH"

echo "[3/3] verifying"
bash "${HERE}/verify_collapse.sh" "$MANIFEST" "$OUTDIR"

echo "wrote ${OUTDIR}/  ($(ls "$OUTDIR"/*.fa 2>/dev/null | wc -l) proteomes)"
