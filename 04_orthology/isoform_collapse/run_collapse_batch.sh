#!/usr/bin/env bash
#
# run_collapse_batch.sh
#
# Reduce a set of NCBI proteomes to one protein per gene, ready for OrthoFinder.
#
# Manifest (TAB-separated, one line per genome, '#' comments allowed):
#     accession <TAB> output_tag <TAB> protein_fasta <TAB> genomic_gff [<TAB> expected_genes]
#
# Build one with make_manifest.py. Two things about the format matter:
#   * fields are separated by TABs, and paths may legitimately contain spaces,
#     so do not align the columns with spaces
#   * save it with Unix line endings; a CR at the end of a line becomes part of
#     expected_genes and argparse rejects it
# The last line may or may not end with a newline; both are handled.
#
# Output is written as <outdir>/<output_tag>.fa -- one extension, so OrthoFinder
# names the species <output_tag> rather than <output_tag>.longest. Any *.fa
# already in <outdir> is removed first, so a shortened manifest cannot leave a
# stale proteome behind to be picked up by the next OrthoFinder run.
#
# Usage:
#     bash run_collapse_batch.sh manifest.tsv clean_proteomes/ [--min-length N]
#
# Then:
#     bash verify_collapse.sh manifest.tsv clean_proteomes/
#     orthofinder -f clean_proteomes/ -t 16
#
set -uo pipefail

MANIFEST="${1:?usage: run_collapse_batch.sh <manifest.tsv> <outdir> [--min-length N]}"
OUTDIR="${2:?usage: run_collapse_batch.sh <manifest.tsv> <outdir> [--min-length N]}"
shift 2
EXTRA=("$@")            # e.g. --min-length 30, passed through to the collapser

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

[[ -f "$MANIFEST" ]] || { echo "no such manifest: $MANIFEST" >&2; exit 2; }

mkdir -p "$OUTDIR"

# Clear previous outputs so a shortened manifest cannot leave a stale proteome
# behind for the next OrthoFinder run to pick up -- but keep the ones this run
# is not going to regenerate, namely genomes the manifest marks pre-collapsed
# with `-` in the genomic_gff field.
keep='|'
while IFS=$'\t' read -r acc tag fasta gff expected || [[ -n "${acc:-}" ]]; do
    tag="${tag%$'\r'}"; gff="${gff%$'\r'}"; acc="${acc%$'\r'}"
    [[ -z "${acc:-}" || "$acc" == \#* ]] && continue
    [[ "$gff" == '-' ]] && keep="${keep}${tag}|"
done < "$MANIFEST"

for f in "$OUTDIR"/*.fa; do
    [[ -e "$f" ]] || continue
    b=$(basename "$f" .fa)
    case "$keep" in
        *"|$b|"*) ;;
        *) rm -f "$f" ;;
    esac
done

fail=0
kept=0
n=0

# `|| [[ -n "$acc" ]]` so a manifest whose last line lacks a trailing newline
# still has that line processed instead of silently dropped.
while IFS=$'\t' read -r acc tag fasta gff expected || [[ -n "${acc:-}" ]]; do
    acc="${acc%$'\r'}"; tag="${tag%$'\r'}"
    fasta="${fasta%$'\r'}"; gff="${gff%$'\r'}"; expected="${expected%$'\r'}"

    [[ -z "${acc:-}" || "$acc" == \#* ]] && continue
    n=$((n + 1))
    echo "============================================================"
    echo "$acc  ->  $tag"

    # `-` in the genomic_gff field marks a genome supplied already collapsed --
    # usually your own assembly. There is nothing to redo, but the file has to
    # be present, because OrthoFinder will read it along with everything else.
    if [[ "$gff" == '-' ]]; then
        if [[ -f "$OUTDIR/${tag}.fa" ]]; then
            echo "  pre-collapsed, left as it is: $OUTDIR/${tag}.fa" \
                 "($(grep -c '^>' "$OUTDIR/${tag}.fa") sequences)"
            kept=$((kept + 1))
        else
            echo "  MISSING: $OUTDIR/${tag}.fa -- the manifest marks this genome" \
                 "pre-collapsed, so the file must already be in the output" \
                 "directory." >&2
            fail=$((fail + 1))
        fi
        continue
    fi

    ok=1
    for f in "$fasta" "$gff"; do
        [[ -f "$f" ]] || { echo "  MISSING: $f" >&2; ok=0; }
    done
    if [[ "$ok" -eq 0 ]]; then
        fail=$((fail + 1))
        continue
    fi

    python3 "$HERE/collapse_isoforms_from_gff.py" \
        --fasta "$fasta" --gff "$gff" \
        --out "$OUTDIR/${tag}.fa" --tag "$tag" \
        ${expected:+--expected-genes "$expected"} \
        ${EXTRA[@]+"${EXTRA[@]}"} || fail=$((fail + 1))
        # ${EXTRA[@]+...} rather than "${EXTRA[@]}": under `set -u`, bash 3.2 --
        # which is what macOS ships -- treats an empty array as unbound.
done < "$MANIFEST"

echo
echo "============================================================"
printf '%-46s %10s\n' proteome sequences
total=0
for f in "$OUTDIR"/*.fa; do
    [[ -e "$f" ]] || continue
    c=$(grep -c '^>' "$f")
    total=$((total + c))
    printf '%-46s %10s\n' "$(basename "$f" .fa)" "$c"
done
printf '%-46s %10s\n' TOTAL "$total"
echo
echo "$n manifest rows, $fail failed, $kept pre-collapsed and left alone"
echo "Now run: bash $(basename "${BASH_SOURCE[0]}" | sed 's/run_collapse_batch/verify_collapse/') \"$MANIFEST\" \"$OUTDIR\""
[[ "$fail" -eq 0 ]] || exit 1
