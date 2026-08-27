#!/usr/bin/env bash
#
# Identification and localization of telomeric repeats with TIDK, used to infer
# chromosome-scale scaffolds. The window table written by `search` is the input
# to 07_telomere_summary.py.
#
# Run from the repository root:
#   conda activate oimp_02_assembly
#   bash 02_assembly/06_telomere_repeats.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
GENOME="${GENOME:-results/02_assembly/fixed/Hifiasm_purged.sorted.fixed.renumbered.fasta}"
OUTDIR="${OUTDIR:-results/02_assembly/tidk}"   # output directory
PREFIX="${PREFIX:-Oimp}"                       # output prefix
MOTIF="${MOTIF:-AACCT}"                        # motif for `search` (TTAGG/AACCT)
CLADE="${CLADE:-Coleoptera}"                   # clade database for `find`
MIN_LEN="${MIN_LEN:-5}"                        # `explore` shortest repeat unit
MAX_LEN="${MAX_LEN:-12}"                       # `explore` longest repeat unit
THRESHOLD="${THRESHOLD:-10}"                   # `explore` minimum repeat run
# -----------------------------------------------------------------------------

for t in tidk seqkit; do
  command -v "$t" >/dev/null || { echo "ERROR: $t not found" >&2; exit 1; }
done
[[ -f "$GENOME" ]] || { echo "ERROR: no such file: $GENOME" >&2; exit 1; }

GENOME=$(readlink -f "$GENOME")
mkdir -p "$OUTDIR"; cd "$OUTDIR"

# 1) Scaffold lengths, needed by 07_telomere_summary.py.
echo "[1/5] scaffold lengths"
seqkit fx2tab -nl "$GENOME" | cut -f1,4 > scaffold_lengths.tsv

# 2) Which simple repeats sit near the scaffold ends? This makes no assumption
#    about the motif. For this assembly the top hits were AAAAAAT, AAAATT and
#    AAAAAT, with the canonical insect repeat AACCT also recovered.
echo "[2/5] tidk explore"
tidk explore --minimum "$MIN_LEN" --maximum "$MAX_LEN" \
             --threshold "$THRESHOLD" --log "$GENOME" \
  > "${PREFIX}_explore_out.tsv"

# 3) Curated Coleoptera motifs (AACCC, AACCT, AACAGACCCG, ACCTG), counted in
#    10 kb windows along every scaffold.
echo "[3/5] tidk find"
tidk find -c "$CLADE" -d . -o "${PREFIX}_find" --log "$GENOME"

# 4) The canonical insect repeat alone, counted the same way. This is the
#    table the telomere calls are made from.
echo "[4/5] tidk search (${MOTIF})"
tidk search --string "$MOTIF" -d . -o "${PREFIX}_search_${MOTIF}" "$GENOME"

# 5) Per-window plots for both.
echo "[5/5] tidk plot"
tidk plot -t "${PREFIX}_find_telomeric_repeat_windows.tsv"            -o "${PREFIX}_find_plot"
tidk plot -t "${PREFIX}_search_${MOTIF}_telomeric_repeat_windows.tsv" -o "${PREFIX}_search_plot_${MOTIF}"

echo "wrote ${OUTDIR}/${PREFIX}_search_${MOTIF}_telomeric_repeat_windows.tsv"
echo "      ${OUTDIR}/scaffold_lengths.tsv"
