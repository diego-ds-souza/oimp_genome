#!/usr/bin/env bash
#
# Build the study gene sets for GO enrichment. Six sets are produced, each a
# list of O. impluviata genes:
#
#   species_specific            orthogroups found only in O. impluviata
#   cerambycidae_specific       orthogroups where every species present is a
#                               cerambycid
#   expanded_vs_Aglabripennis   more genes than in A. glabripennis
#   expanded_vs_Cerambycidae    more genes than the maximum of the other
#                               Cerambycidae, O. impluviata excluded
#   expanded_vs_Chrysomelidae   more genes than the maximum Chrysomelidae
#   expanded_vs_Curculionoidea  more genes than the maximum Curculionoidea
#
# The reference for every expansion set is the MAXIMUM of the group, not its
# mean, so a set counts only orthogroups larger in O. impluviata than in any
# member of the comparison group. Clade membership comes from
# species_taxonomy.tsv, written by 04_orthology/07_root_tree.R.
#
# Run from the repository root:
#   conda activate oimp_06_go_enrichment
#   bash 06_go_enrichment/00_build_gene_sets.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
PROJECT="${PROJECT:-results/04_orthology}"       # where OrthoFinder wrote
OF_BASE="${OF_BASE:-${PROJECT}/orthofinder_out}"
TAXONOMY="${TAXONOMY:-${PROJECT}/species_taxonomy.tsv}"
OUT="${OUT:-results/06_go_enrichment}"           # output directory
OI_KEY="${OI_KEY:-Oncideres_impluviata}"         # focal species column
AG_KEY="${AG_KEY:-Anoplophora_glabripennis}"     # single-species reference
MIN_FOLD_CHANGE="${MIN_FOLD_CHANGE:-1.0}"        # 1.0 = any increase counts
MIN_GENE_DIFF="${MIN_GENE_DIFF:-1}"              # minimum extra genes
# -----------------------------------------------------------------------------

# 1) Locate the OrthoFinder results, preferring the path step 03 recorded.
if [[ -f "${OF_BASE}/LATEST_RESULTS_PATH.txt" ]]; then
  RES=$(cat "${OF_BASE}/LATEST_RESULTS_PATH.txt")
else
  RES=$(ls -d ${OF_BASE}/run_*/Results_* 2>/dev/null | sort | tail -n1)
fi

OG_TSV="${RES}/Orthogroups/Orthogroups.tsv"
GC_TSV="${RES}/Orthogroups/Orthogroups.GeneCount.tsv"

mkdir -p "$OUT"

echo "[CONFIG] MIN_FOLD_CHANGE = ${MIN_FOLD_CHANGE}"
echo "[CONFIG] MIN_GENE_DIFF = ${MIN_GENE_DIFF}"

python3 - << 'PY' "$OG_TSV" "$GC_TSV" "$OUT" "$OI_KEY" "$AG_KEY" "$MIN_FOLD_CHANGE" "$MIN_GENE_DIFF" "$TAXONOMY"
import sys, pandas as pd, re

og_tsv, gc_tsv, out, OI_KEY, AG_KEY = sys.argv[1:6]
MIN_FOLD_CHANGE = float(sys.argv[6])
MIN_GENE_DIFF = int(sys.argv[7])
taxonomy_file = sys.argv[8]

# --- Load the two Orthofinder tables
og = pd.read_csv(og_tsv, sep="\t", dtype=str).fillna("")
gc = pd.read_csv(gc_tsv, sep="\t")
species_counts_cols = [c for c in gc.columns if c not in ("Orthogroup","Total")]

# Resolve exact column names. Species names are now plain genus_species, so an
# exact match is expected; the containment fallback is kept for safety.
def resolve_col(key, cols):
    if key in cols:
        return key
    cands = [c for c in cols if key in c]
    if len(cands)==1:
        return cands[0]
    if len(cands)>1:
        print(f"[WARNING] Multiple matches for {key}: {cands}")
        return cands[0]
    return None

OI = resolve_col(OI_KEY, list(og.columns))
AG = resolve_col(AG_KEY, species_counts_cols)

if OI is None:
    raise SystemExit(f"ERROR: Could not resolve O. impluviata column using key '{OI_KEY}'. Available columns include: {list(og.columns[:10])} ...")
if AG is None:
    raise SystemExit(f"ERROR: Could not resolve A. glabripennis column using key '{AG_KEY}'. Available species columns include: {species_counts_cols[:10]} ...")

print(f"[INFO] Using O. impluviata column: {OI}")
print(f"[INFO] Using A. glabripennis column: {AG}")

# Keep counts table indexed by OG
gc = gc[gc["Orthogroup"]!="Total"].copy().set_index("Orthogroup")
species = species_counts_cols

# -----------------------------------------------------------------------------
# Clade membership from species_taxonomy.tsv
#
# Replaces the v4 filename-prefix parsing (CER_ / CHRY_ / CUR_ / BREN_ / ...).
# The assignment rules are the same: Curculionoidea is checked at superfamily
# level first, then Cerambycidae and Chrysomelidae at family level. Anything
# else (outgroups, unplaced taxa) falls through to "Other" and is excluded from
# every reference group, exactly as "Outgroup"/"Unknown" were in v4.
# -----------------------------------------------------------------------------
taxonomy = pd.read_csv(taxonomy_file, sep="\t", dtype=str).fillna("")
required = {"Species","Family","Superfamily"}
missing_cols = required - set(taxonomy.columns)
if missing_cols:
    raise SystemExit("ERROR: missing taxonomy columns: " + ", ".join(sorted(missing_cols)))
taxonomy = taxonomy.set_index("Species")

def clade_of(sp):
    if sp not in taxonomy.index:
        return "Unknown"
    family = taxonomy.loc[sp, "Family"]
    superfamily = taxonomy.loc[sp, "Superfamily"]
    if superfamily == "Curculionoidea": return "Curculionoidea"
    if family == "Cerambycidae":        return "Cerambycidae"
    if family == "Chrysomelidae":       return "Chrysomelidae"
    return "Other"

clades = {sp:clade_of(sp) for sp in species}

unknown = [sp for sp in species if clades[sp]=="Unknown"]
if unknown:
    print(f"[WARNING] {len(unknown)} species missing from {taxonomy_file}: {unknown}")
    print("[WARNING] They are excluded from all reference groups.")

curc = [s for s in species if clades[s]=="Curculionoidea"]
chry = [s for s in species if clades[s]=="Chrysomelidae"]
cer  = [s for s in species if clades[s]=="Cerambycidae"]

# Other Cerambycidae (excluding O. impluviata itself)
OI_col_in_counts = resolve_col(OI_KEY, species_counts_cols)
if OI_col_in_counts is None:
    raise SystemExit(f"ERROR: Could not resolve O. impluviata in the GeneCount columns using key '{OI_KEY}'")
cer_other = [s for s in cer if s != OI_col_in_counts]

print(f"[INFO] Cerambycidae species (n={len(cer)}): {cer}")
print(f"[INFO] Other Cerambycidae for comparison (n={len(cer_other)}): {cer_other}")
print(f"[INFO] Chrysomelidae species (n={len(chry)})")
print(f"[INFO] Curculionoidea species (n={len(curc)})")

OGs = list(gc.index)
OGset = set(OGs)

def present_species(ogid):
    row = gc.loc[ogid]
    return {sp for sp in species if int(row[sp])>0}

# Map OG -> O_impluviata genes using Orthogroups.tsv column
OG2GENES_OI = {}
for _, row in og.iterrows():
    ogid = row["Orthogroup"]
    if ogid not in OGset:
        continue
    li = row[OI].strip()
    if not li:
        continue
    genes = [g.strip() for g in re.split(r'[,\s]+', li) if g.strip()]
    if genes:
        OG2GENES_OI[ogid] = genes

# -----------------------------------------------------------------------------
# Helper function for expansion criteria
# -----------------------------------------------------------------------------
def is_expanded(oi_count, ref_count, min_fold=MIN_FOLD_CHANGE, min_diff=MIN_GENE_DIFF):
    """
    Check if O. impluviata gene count represents an expansion relative to reference.

    Args:
        oi_count: Gene count in O. impluviata
        ref_count: Gene count in reference (single species or max of group)
        min_fold: Minimum fold-change required (default from config)
        min_diff: Minimum absolute difference required (default from config)

    Returns:
        True if expanded, False otherwise
    """
    if oi_count <= ref_count:
        return False

    # Check absolute difference
    if (oi_count - ref_count) < min_diff:
        return False

    # Check fold-change (handle ref_count == 0)
    if ref_count == 0:
        # If reference has 0 and OI has genes, it's an expansion
        return oi_count >= min_diff

    fold = oi_count / ref_count
    return fold >= min_fold

# -----------------------------------------------------------------------------
# Study sets
# -----------------------------------------------------------------------------

# A. Species-specific: OGs found only in O. impluviata
A_species_specific = [ogid for ogid in OGs if present_species(ogid)=={OI_col_in_counts}]

# B. Cerambycidae-specific: OGs found only in Cerambycidae
B_cer_only = [ogid for ogid in OGs
              if all(sp in cer for sp in present_species(ogid))
              and not any(sp in chry for sp in present_species(ogid))
              and not any(sp in curc for sp in present_species(ogid))]

# C. Expanded vs A. glabripennis (single species comparison)
C_exp_vs_AG = []
for ogid in OGs:
    oi = int(gc.loc[ogid, OI_col_in_counts])
    ag = int(gc.loc[ogid, AG])
    if is_expanded(oi, ag):
        C_exp_vs_AG.append(ogid)

# D. Expanded vs other Cerambycidae (max of other cerambycids)
D_exp_vs_CER = []
if cer_other:
    for ogid in OGs:
        oi = int(gc.loc[ogid, OI_col_in_counts])
        mx = max(int(gc.loc[ogid, s]) for s in cer_other)
        if is_expanded(oi, mx):
            D_exp_vs_CER.append(ogid)
else:
    print("[WARN] No other Cerambycidae species found for vs_Cerambycidae comparison")

# E. Expanded vs Chrysomelidae (max of all chrysomelids)
E_exp_vs_CHRY = []
if chry:
    for ogid in OGs:
        oi = int(gc.loc[ogid, OI_col_in_counts])
        mx = max(int(gc.loc[ogid, s]) for s in chry)
        if is_expanded(oi, mx):
            E_exp_vs_CHRY.append(ogid)

# F. Expanded vs Curculionoidea (max of all curculionoids)
F_exp_vs_CURC = []
if curc:
    for ogid in OGs:
        oi = int(gc.loc[ogid, OI_col_in_counts])
        mx = max(int(gc.loc[ogid, s]) for s in curc)
        if is_expanded(oi, mx):
            F_exp_vs_CURC.append(ogid)

# -----------------------------------------------------------------------------
# Write gene lists
# -----------------------------------------------------------------------------
def clean_id(g):
    """Strip species prefix and isoform suffix from gene ID."""
    core = g.split("|")[-1]
    core = core.split("-")[0]
    return core

def write_set(name, og_list):
    """Write study set to file and report statistics."""
    genes=[]
    for ogid in og_list:
        genes += OG2GENES_OI.get(ogid, [])
    uniq = sorted({clean_id(x) for x in genes})
    with open(f"{out}/study_{name}.txt","w") as h:
        for x in uniq: h.write(x+"\n")
    print(f"{name:30s} OGs={len(og_list):6d}  genes={len(uniq):6d}")

print("\n" + "="*70)
print("STUDY SETS SUMMARY")
print("="*70)

SETS = [
    ("species_specific",           A_species_specific),
    ("cerambycidae_specific",      B_cer_only),
    ("expanded_vs_Aglabripennis",  C_exp_vs_AG),
    ("expanded_vs_Cerambycidae",   D_exp_vs_CER),
    ("expanded_vs_Chrysomelidae",  E_exp_vs_CHRY),
    ("expanded_vs_Curculionoidea", F_exp_vs_CURC),
]

for nm, ogs in SETS:
    write_set(nm, ogs)

# -----------------------------------------------------------------------------
# Write summary statistics file for methods reporting
# -----------------------------------------------------------------------------
summary_file = f"{out}/study_sets_summary.tsv"
with open(summary_file, "w") as sf:
    sf.write("study_set\tn_orthogroups\tn_genes\tmin_fold_change\tmin_gene_diff\n")
    for nm, ogs in SETS:
        genes = []
        for ogid in ogs:
            genes += OG2GENES_OI.get(ogid, [])
        n_genes = len({clean_id(x) for x in genes})
        sf.write(f"{nm}\t{len(ogs)}\t{n_genes}\t{MIN_FOLD_CHANGE}\t{MIN_GENE_DIFF}\n")

print(f"\nSummary written to: {summary_file}")
PY

echo ""
echo "Study gene sets written to ${OUT}/study_*.txt"
echo ""
echo "To use more stringent expansion criteria, run with:"
echo "  MIN_FOLD_CHANGE=1.5 MIN_GENE_DIFF=2 bash $0"
