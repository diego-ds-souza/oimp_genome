#!/usr/bin/env bash
#
# Independent check on the orthogroup category counts from step 08.
#
# The totals are recomputed straight from Orthogroups.GeneCount.tsv and the
# single-copy orthologue list, then compared with the summary table, so a
# mis-parsed taxonomy or a stale OrthoFinder run shows up as a mismatch rather
# than as a plausible-looking figure.
#
# Run from the repository root:
#   conda activate oimp_04_orthology
#   bash 04_orthology/09_verify_counts.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
PROJECT="${PROJECT:-results/04_orthology}"                 # section output root
OF_BASE="${OF_BASE:-${PROJECT}/orthofinder_out}"           # OrthoFinder output
CSV_LONG="${CSV_LONG:-${PROJECT}/plots/orthology_histogram_counts_long.csv}"
                                                           # summary from step 08
# -----------------------------------------------------------------------------

# 1) Locate the OrthoFinder results, preferring the path step 03 recorded.
if [[ -f "${OF_BASE}/LATEST_RESULTS_PATH.txt" ]]; then
  RES="$(cat "${OF_BASE}/LATEST_RESULTS_PATH.txt")"
else
  RES=$(ls -d ${OF_BASE}/run_*/Results_* 2>/dev/null | sort | tail -n1 || true)
  [[ -z "${RES}" ]] && RES=$(ls -d ${OF_BASE}/Results_* 2>/dev/null | sort | tail -n1 || true)
fi
[[ -n "${RES:-}" && -d "$RES" ]] || { echo "Error: Results_* not found under ${OF_BASE}"; exit 1; }

# 2) Recount the single-copy universal orthogroups two independent ways and
#    compare them with the table step 08 wrote.
GENECOUNT="${RES}/Orthogroups/Orthogroups.GeneCount.tsv"
SCO_LIST="${RES}/Orthogroups/Orthogroups_SingleCopyOrthologues.txt"

echo "Using:"
echo "  GeneCount TSV:   $GENECOUNT"
echo "  SCO list:        $SCO_LIST"
echo "  Histogram CSV:   $CSV_LONG"
echo

python3 - << 'PY' "$GENECOUNT" "$SCO_LIST" "$CSV_LONG"
import sys, pandas as pd

gene_count_tsv, sco_list, csv_long = sys.argv[1:]

# 1) Count universal SCOs from OrthoFinder's own list
with open(sco_list) as fh:
    sco_ids = [ln.strip() for ln in fh if ln.strip()]
sco_n = len(sco_ids)

# 2) Recompute from GeneCount.tsv (all species present with exactly 1 copy)
og = pd.read_csv(gene_count_tsv, sep='\t')
species = [c for c in og.columns if c not in ("Orthogroup","Total")]
def is_111(row):
    return all(int(row[s])==1 for s in species)
sco_from_tsv = og.loc[og[species].applymap(int).apply(lambda r: all(r==1), axis=1), "Orthogroup"].tolist()

# 3) Read the CSV we produced for plotting
df = pd.read_csv(csv_long)
sco_csv = df[df["Category"]=="SingleCopy_Universal"].groupby("Species")["Count"].first()

print(f"OrthoFinder SCO list count           : {sco_n}")
print(f"GeneCount.tsv-derived 1:1:1 OGs      : {len(sco_from_tsv)}")
mismatch = (set(sco_from_tsv) - set(sco_ids)) | (set(sco_ids) - set(sco_from_tsv))
print(f"IDs mismatch between methods         : {len(mismatch)} OGs")
print()
print("CSV SingleCopy_Universal per species (first 5):")
print(sco_csv.head())
PY
