#!/usr/bin/env bash
#
# verify_collapse.sh -- the check that catches a silently ineffective collapse.
#
# Usage: bash verify_collapse.sh manifest.tsv clean_dir
#
# It reads the raw proteome paths and the expected gene counts straight out of
# the manifest, so it works on an ncbi-datasets download as it comes off the
# wire -- where every proteome is called protein.faa and sits in its own
# accession folder, and no directory of per-species raw files exists.
#
# The primary test is arithmetic: the collapsed sequence count must equal the
# number of protein-coding genes in the annotation those proteins came from.
# That is exact, and it holds for RefSeq and GenBank alike.
#
# It also reports any proteome sitting in the clean directory that the manifest
# does not mention. OrthoFinder reads the whole directory, so such a file is
# included in the analysis while being covered by no check at all.
#
# The ` isoform X` column is reported for information only. It is a RefSeq
# naming convention; GenBank submissions carry redundant transcripts with no
# such marker, so a check built on it passes a failed collapse on GenBank
# assemblies without comment. Treat the column as a description of the input,
# never as the test.
#
# Exit status is non-zero if any proteome is flagged.
#
set -uo pipefail

MANIFEST="${1:?usage: verify_collapse.sh <manifest.tsv> <clean_dir>}"
CLEAN="${2:?usage: verify_collapse.sh <manifest.tsv> <clean_dir>}"

[[ -f "$MANIFEST" ]] || { echo "no such manifest: $MANIFEST" >&2; exit 2; }
[[ -d "$CLEAN"    ]] || { echo "no such directory: $CLEAN" >&2; exit 2; }

printf '%-32s %9s %9s %9s %9s %11s  %s\n' \
       proteome raw clean removed expected 'isoform X' status

flagged=0
unchecked=0
asserted=0
rows=0
seen='|'          # tags covered by the manifest, for the orphan check below

while IFS=$'\t' read -r acc tag fasta gff expected || [[ -n "${acc:-}" ]]; do
    acc="${acc%$'\r'}"; tag="${tag%$'\r'}"
    fasta="${fasta%$'\r'}"; gff="${gff%$'\r'}"; expected="${expected%$'\r'}"
    [[ -z "${acc:-}" || "$acc" == \#* ]] && continue
    rows=$((rows + 1))
    seen="${seen}${tag}|"

    # A genome supplied already collapsed -- typically your own assembly, whose
    # annotation this pipeline never processed -- is written with `-` in the
    # genomic_gff field. Its expected_genes is then an assertion by whoever
    # wrote the manifest rather than a count taken from an annotation, so it is
    # checked against the file but labelled to say where the number came from.
    precollapsed=0
    [[ "$gff" == '-' ]] && precollapsed=1

    if [[ "$precollapsed" -eq 0 && ! -f "$fasta" ]]; then
        printf '%-32s %9s  %s\n' "$tag" - "*** RAW FASTA MISSING ***"
        flagged=$((flagged + 1)); continue
    fi

    c=""
    for cand in "$CLEAN/$tag.fa" "$CLEAN/$tag.faa" "$CLEAN/$tag.longest.faa"; do
        [[ -f "$cand" ]] && { c="$cand"; break; }
    done
    if [[ -z "$c" ]]; then
        printf '%-32s %9s  %s\n' "$tag" - "*** NO CLEAN FILE ***"
        flagged=$((flagged + 1)); continue
    fi

    nc=$(grep -c '^>' "$c" || true)
    if [[ -f "$fasta" ]]; then
        nr=$(grep -c '^>' "$fasta" || true)
        iso=$(grep -c ' isoform X' "$fasta" || true)
    else
        nr=""; iso=""      # pre-collapsed row with no source proteome recorded
    fi

    if [[ "$precollapsed" -eq 1 ]]; then
        if [[ -z "$expected" ]]; then
            st="*** PRE-COLLAPSED, NO expected_genes ***"; flagged=$((flagged + 1))
            expected="-"
        elif [[ "$nc" -eq "$expected" ]]; then
            st="ok (pre-collapsed; gene count asserted, not derived)"
            asserted=$((asserted + 1))
        else
            st="*** CHECK ***  $((nc - expected)) vs asserted"
            flagged=$((flagged + 1))
        fi
        printf '%-32s %9s %9s %9s %9s %11s  %s\n' \
               "$tag" "${nr:--}" "$nc" "-" "$expected" "${iso:--}" "$st"
        continue
    fi

    if [[ -z "$expected" ]]; then
        # No expected count: fall back to the weak heuristic and say so.
        if [[ "$nr" -eq "$nc" && "$iso" -gt 0 ]]; then
            st="*** NOTHING COLLAPSED ***"; flagged=$((flagged + 1))
        else
            st="unchecked (no expected_genes)"; unchecked=$((unchecked + 1))
        fi
        expected="-"
    elif [[ "$nc" -eq "$expected" ]]; then
        st="ok"
    elif [[ "$nr" -eq "$nc" ]]; then
        st="*** NOTHING COLLAPSED ***"; flagged=$((flagged + 1))
    else
        st="*** CHECK ***  $((nc - expected)) vs expected"
        flagged=$((flagged + 1))
    fi

    printf '%-32s %9s %9s %9s %9s %11s  %s\n' \
           "$tag" "$nr" "$nc" "$((nr - nc))" "$expected" "$iso" "$st"
done < "$MANIFEST"

# Anything sitting in the clean directory that the manifest does not mention.
# OrthoFinder will happily include such a file -- it reads the whole directory --
# and nothing above will have looked at it, so it is exactly the kind of
# unchecked proteome this script exists to make impossible. Usually it is a
# genome added by hand: your own assembly, or one from outside NCBI.
orphans=0
for f in "$CLEAN"/*.fa "$CLEAN"/*.faa; do
    [[ -e "$f" ]] || continue
    b=$(basename "$f"); b="${b%.gz}"; b="${b%.fa}"; b="${b%.faa}"; b="${b%.longest}"
    case "$seen" in
        *"|$b|"*) ;;
        *)
            if [[ "$orphans" -eq 0 ]]; then
                echo
                echo "NOT IN THE MANIFEST -- present in $CLEAN but never checked:"
            fi
            n=$(grep -c '^>' "$f" || true)
            printf '  %-32s %9s sequences\n' "$b" "$n"
            orphans=$((orphans + 1))
            ;;
    esac
done
if [[ "$orphans" -gt 0 ]]; then
    echo "  OrthoFinder reads the whole directory, so these WILL be included."
    echo "  Add each one to the manifest -- with its own annotation GFF and"
    echo "  protein-coding gene count -- so it gets the same check as the rest."
fi

echo
checked=$((rows - unchecked))
if [[ "$flagged" -gt 0 || "$orphans" -gt 0 ]]; then
    echo "$rows proteomes from the manifest, $flagged FLAGGED," \
         "$orphans not in the manifest." >&2
    exit 1
elif [[ "$unchecked" -gt 0 ]]; then
    echo "$rows proteomes: $checked matched their protein-coding gene count," \
         "$unchecked NOT CHECKED (no expected_genes in the manifest)." >&2
    echo "Rebuild the manifest with make_manifest.py to get a real check." >&2
    exit 1
elif [[ "$asserted" -gt 0 ]]; then
    echo "$rows proteomes checked, all matched their expected gene count" \
         "-- $((rows - asserted)) derived from an annotation," \
         "$asserted asserted in the manifest (pre-collapsed)."
else
    echo "$rows proteomes checked, all matched their protein-coding gene count."
fi
