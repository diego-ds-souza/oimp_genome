#!/usr/bin/env python3
"""
GO enrichment for each study set, with GOAtools.

Fisher's exact test, counts propagated up the GO hierarchy
(propagate_counts = True), and Benjamini-Hochberg FDR correction. Terms with
FDR < 0.05 are reported as significant. Each namespace - biological process,
molecular function, cellular component - is tested separately, as is each of
the six comparative baselines.

Run from the repository root:
    conda activate oimp_06_go_enrichment
    python 06_go_enrichment/02_go_enrichment.py
"""
import os
import sys
import pandas as pd
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy

# --------------------------------- settings ----------------------------------
SETS_DIR = os.environ.get("SETS_DIR", "results/06_go_enrichment")
OBO = os.environ.get("OBO", os.path.join(SETS_DIR, "go-basic.obo"))
GENE2GO = os.environ.get("GENE2GO", os.path.join(SETS_DIR, "Oi_gene2go.tsv"))
OUT_DIR = os.environ.get("OUT_DIR", os.path.join(SETS_DIR, "results"))
ALPHA = float(os.environ.get("ALPHA", 0.05))     # FDR threshold
# -----------------------------------------------------------------------------

os.makedirs(OUT_DIR, exist_ok=True)

# The six comparative baselines, from step 00.
studies = [
    ("species_specific",           f"{SETS_DIR}/study_species_specific.txt"),
    ("cerambycidae_specific",      f"{SETS_DIR}/study_cerambycidae_specific.txt"),
    ("expanded_vs_Aglabripennis",  f"{SETS_DIR}/study_expanded_vs_Aglabripennis.txt"),
    ("expanded_vs_Cerambycidae",   f"{SETS_DIR}/study_expanded_vs_Cerambycidae.txt"),
    ("expanded_vs_Chrysomelidae",  f"{SETS_DIR}/study_expanded_vs_Chrysomelidae.txt"),
    ("expanded_vs_Curculionoidea", f"{SETS_DIR}/study_expanded_vs_Curculionoidea.txt"),
]

print("="*70)
print("GO Enrichment Analysis (GOAtools)")
print("="*70)

# Load GO DAG
print(f"\nLoading GO ontology: {OBO}")
obodag = GODag(OBO, optional_attrs={'relationship'})

# Load gene2go mapping
print(f"Loading gene-GO associations: {GENE2GO}")
gene2gos = {}
with open(GENE2GO) as f:
    for ln in f:
        if not ln.strip():
            continue
        gid, go = ln.strip().split("\t")
        gene2gos.setdefault(gid, set()).add(go)

pop_ids = set(gene2gos.keys())
print(f"Background population: {len(pop_ids)} genes with GO annotations")

# Initialize GOAtools enrichment object
goeaobj = GOEnrichmentStudy(
    pop_ids,
    gene2gos,
    obodag,
    propagate_counts=True,
    alpha=ALPHA,
    methods=['fdr_bh']
)

print("\n" + "="*70)
print("Running enrichment tests")
print("="*70)

# Summary for all studies
summary_rows = []

for label, path in studies:
    print(f"\n--- {label} ---")
    
    if not os.path.isfile(path):
        print(f"[WARN] Missing study set: {path}; skipping")
        continue
    
    with open(path) as f:
        study_ids_raw = [x.strip() for x in f if x.strip()]
    
    # Intersect with background (only genes with GO annotations)
    study_ids = [g for g in study_ids_raw if g in pop_ids]
    n_no_go = len(study_ids_raw) - len(study_ids)
    
    if not study_ids:
        print(f"[WARN] No study genes with GO annotations; skipping")
        continue
    
    print(f"  Study genes: {len(study_ids_raw)} total, {len(study_ids)} with GO ({n_no_go} without)")
    
    # Run enrichment
    res = goeaobj.run_study(study_ids)
    
    # Build results table with additional metrics
    rows = []
    for r in res:
        # Calculate fold-enrichment
        study_ratio = r.study_count / len(study_ids) if len(study_ids) > 0 else 0
        pop_ratio = r.pop_count / len(pop_ids) if len(pop_ids) > 0 else 0
        fold_enrichment = study_ratio / pop_ratio if pop_ratio > 0 else float('inf')
        
        rows.append({
            "GO": r.GO,
            "Term": {"BP": "BP", "MF": "MF", "CC": "CC"}[r.NS],
            "Name": r.name,
            "study_count": r.study_count,
            "study_total": len(study_ids),
            "pop_count": r.pop_count,
            "pop_total": len(pop_ids),
            "Ratio in study": f"{r.study_count}/{len(study_ids)}",
            "Ratio in pop": f"{r.pop_count}/{len(pop_ids)}",
            "fold_enrichment": round(fold_enrichment, 3),
            "pvalue": r.p_uncorrected,
            "p.adjust": r.p_fdr_bh
        })
    
    df = pd.DataFrame(rows)
    
    # Sort by adjusted p-value, then raw p-value
    df = df.sort_values(["p.adjust", "pvalue"])
    
    # Save FULL results (for supplementary materials)
    out_full = f"{OUT_DIR}/GO_{label}_full.tsv"
    df.to_csv(out_full, sep="\t", index=False)
    
    # Save FILTERED results (significant only, for main analysis)
    df_sig = df[df["p.adjust"] < ALPHA].copy()
    out_sig = f"{OUT_DIR}/GO_{label}.tsv"
    df_sig.to_csv(out_sig, sep="\t", index=False)
    
    # Count significant terms by namespace
    n_sig_bp = len(df_sig[df_sig["Term"] == "BP"])
    n_sig_mf = len(df_sig[df_sig["Term"] == "MF"])
    n_sig_cc = len(df_sig[df_sig["Term"] == "CC"])
    n_sig_total = len(df_sig)
    
    print(f"  Significant GO terms (FDR < {ALPHA}): {n_sig_total} (BP={n_sig_bp}, MF={n_sig_mf}, CC={n_sig_cc})")
    print(f"  Output: {out_sig}")
    print(f"  Full results: {out_full}")
    
    # Add to summary
    summary_rows.append({
        "study_set": label,
        "n_study_genes_total": len(study_ids_raw),
        "n_study_genes_with_GO": len(study_ids),
        "n_background_genes": len(pop_ids),
        "n_significant_BP": n_sig_bp,
        "n_significant_MF": n_sig_mf,
        "n_significant_CC": n_sig_cc,
        "n_significant_total": n_sig_total
    })

# Write summary table
summary_df = pd.DataFrame(summary_rows)
summary_file = f"{OUT_DIR}/GO_enrichment_summary.tsv"
summary_df.to_csv(summary_file, sep="\t", index=False)

print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(summary_df.to_string(index=False))
print(f"\nSummary saved to: {summary_file}")
print("\nDone!")
