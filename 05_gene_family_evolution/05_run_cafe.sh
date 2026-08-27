#!/usr/bin/env bash
#
# Runs CAFE5 on the filtered gene family counts and summarises the result.
#
# Four models are fitted, in increasing order of realism:
#   All runs use -p, which estimates a Poisson root-size distribution. Without
#   it CAFE5 assumes a uniform distribution over root sizes, which underflows to
#   -lnL = inf when the data contain families with a large max-min spread.
#
#   base        one birth-death rate (lambda) for the whole tree
#   base_error  the same, with an assembly/annotation error model estimated
#               from the data (-e). This is the run to report: it prevents
#               annotation noise being read as gene family change.
#   gamma_k3    rate heterogeneity across families, three gamma categories
#   large       the extreme families that were split out during preparation
#
# Inputs and outputs both live in $RES, so a run is self-contained and can be
# copied elsewhere in one directory.
#
# The four runs together take several hours; THREADS=20 helps.
#
# Run from the repository root:
#   conda activate oimp_05_gene_family_evolution
#   bash 05_gene_family_evolution/05_run_cafe.sh
#
set -uo pipefail

# ------------------------------- settings ------------------------------------
RES="${RES:-results/05_gene_family_evolution}"   # inputs and outputs both live here
COUNTS="${COUNTS:-${RES}/cafe_gene_counts_filtered.tsv}"
COUNTS_LARGE="${COUNTS_LARGE:-${RES}/cafe_gene_counts_large.tsv}"
TREE="${TREE:-${RES}/species_tree_ultrametric_fixed.nwk}"
THREADS="${THREADS:-16}"                         # CPU threads
# -----------------------------------------------------------------------------

# 1) Check CAFE5 is available and the inputs exist.
command -v cafe5 >/dev/null || {
  echo "ERROR: cafe5 not found. conda activate cafe5" >&2; exit 1; }

for f in "$COUNTS" "$TREE"; do
  [[ -f "$f" ]] || { echo "ERROR: missing input $f" >&2; exit 1; }
done

echo "CAFE5   : $(cafe5 --version 2>&1 | head -1)"
echo "counts  : $COUNTS  ($(( $(wc -l < "$COUNTS") - 1 )) families)"
echo "tree    : $TREE"
echo "outdir  : $RES"
echo "threads : $THREADS"
echo

# 2) The tree tips and the count-table columns must name the same species;
#    CAFE5 matches by name and any mismatch is fatal.
python3 - "$COUNTS" "$TREE" <<'PY' || exit 1
import re, sys
cols = open(sys.argv[1]).readline().rstrip('\n').split('\t')[2:]
nwk  = open(sys.argv[2]).read()
tips = set(re.findall(r'[(,]\s*([A-Za-z0-9_.\-]+)\s*:', nwk))
missing_in_tree  = [c for c in cols if c not in tips]
missing_in_table = [t for t in tips if t not in cols]
if missing_in_tree or missing_in_table:
    print('FATAL: species names do not match.')
    if missing_in_tree:  print('  in table but not in tree:', ', '.join(missing_in_tree))
    if missing_in_table: print('  in tree but not in table:', ', '.join(missing_in_table))
    sys.exit(1)
print(f'species names match ({len(cols)} taxa)')
PY

# 3) The two numerical traps that make every likelihood inf: a family whose
#    max-min spread is too wide, and a branch that is numerically zero.
python3 - "$COUNTS" "$TREE" <<'PY' || exit 1
import csv, re, sys

# 1. family size spread
rows = list(csv.reader(open(sys.argv[1]), delimiter='\t'))
d = sorted(((max(v) - min(v), r[1]) for r in rows[1:]
            for v in [[int(x) for x in r[2:]]]), reverse=True)
print(f'largest max-min spread : {d[0][0]} ({d[0][1]})')
if d[0][0] > 40:
    print('  FATAL: CAFE5 commonly fails to initialise above ~40.')
    print('  Re-run 01_prepare_cafe_input.py with --max-count 30 (or 20).')
    sys.exit(1)

# 2. numerically tiny branches. A branch that is positive but below ~1e-6 is
#    the failure mode that is easy to miss: it is not zero and not negative,
#    so ape's is.ultrametric() and a simple <= 0 test both pass it.
lens = [float(x) for x in re.findall(r':([0-9.eE+-]+)', open(sys.argv[2]).read())]
tiny = [L for L in lens if L < 1e-6]
print(f'shortest branch        : {min(lens):g}')
if tiny:
    print(f'  FATAL: {len(tiny)} branch(es) below 1e-6.')
    print('  Run 02_b_fix_tree.py --min-branch 0.25 and use the _fixed tree.')
    sys.exit(1)
PY
echo

# 4) Fit the four models.
run () {                       # run <outdir> <extra cafe5 args...>
  local out="$1"; shift
  local dir="${RES}/${out}"
  echo "=== ${out}"
  if compgen -G "${dir}/*_family_results.txt" >/dev/null; then
    echo "    already done, skipping (delete ${dir} to force a re-run)"
    return 0
  fi
  cafe5 -i "$COUNTS" -t "$TREE" -c "$THREADS" -o "$dir" "$@" \
      > "${RES}/${out}.log" 2>&1 \
    && echo "    done" \
    || { echo "    FAILED - see ${RES}/${out}.log"; return 1; }
}

run base       -p
run base_error -p -e
run gamma_k3   -p -k 3

echo "=== large"
if [[ -s "$COUNTS_LARGE" && $(wc -l < "$COUNTS_LARGE") -gt 1 ]]; then
  if compgen -G "${RES}/large/*_family_results.txt" >/dev/null; then
    echo "    already done, skipping"
  else
    cafe5 -i "$COUNTS_LARGE" -t "$TREE" -c "$THREADS" -p -o "${RES}/large" \
        > "${RES}/large.log" 2>&1 && echo "    done" \
      || echo "    FAILED - see ${RES}/large.log"
  fi
else
  echo "    no large families, skipped"
fi

# 5) Summarise what each run produced.
echo
echo "================ summary ================"
for d in base base_error gamma_k3 large; do
  res="${RES}/${d}"
  fam=$(ls "$res"/*_family_results.txt 2>/dev/null | head -1)
  [[ -n "$fam" ]] || continue
  echo "${d}:"
  # Lambda / Epsilon / Alpha all live in *_results.txt, one per line
  grep -h -E "^(Model|Lambda|Epsilon|Alpha)" "$res"/*_results.txt 2>/dev/null \
    | sed 's/^/   /'
  total=$(( $(wc -l < "$fam") - 1 ))
  # column 3 is CAFE's own "Significant at 0.05" y/n flag
  sig=$(awk -F'\t' 'NR>1 && $3=="y"' "$fam" | wc -l)
  echo "   families tested            : ${total}"
  echo "   significant (p < 0.05)     : ${sig}"
  clade=$(ls "$res"/*_clade_results.txt 2>/dev/null | head -1)
  [[ -n "$clade" ]] && echo "   per-branch counts          : ${clade}"
done

echo
echo "Report the base_error run. Files in ${RES}/base_error/:"
echo "  Base_results.txt              lambda, epsilon, final -lnL"
echo "  Base_report.cafe              the tree with CAFE's node IDs"
echo "  Base_family_results.txt       one row per family, family-wide p-value"
echo "  Base_clade_results.txt        expansions and contractions per branch"
echo "  Base_change.tab               inferred size change on every branch"
echo "  Base_branch_probabilities.tab which branch drives each significant family"
echo "  Base_asr.tre                  ancestral family sizes, annotated Newick"
