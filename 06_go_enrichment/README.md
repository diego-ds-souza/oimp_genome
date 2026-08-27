# 06 — GO enrichment

Methods section 6, final part: testing whether the orthogroups that expanded in
*O. impluviata* are biased toward particular functions (Fig. 3).

## Scripts

| Script | Software | Purpose |
|---|---|---|
| `00_build_gene_sets.sh` | — | Build the six study gene sets |
| `01_gene2go_from_eggnog.sh` | — | Gene-to-GO background from eggNOG-mapper |
| `02_go_enrichment.py` | GOAtools | Fisher's exact test with FDR correction |
| `03_plot_go_panels.R` | ggplot2 | Dot plots (Fig. 3) and the term tables |

## Running

```bash
conda env create -f 06_go_enrichment/environment.yml
conda activate oimp_06_go_enrichment

bash    06_go_enrichment/00_build_gene_sets.sh
bash    06_go_enrichment/01_gene2go_from_eggnog.sh
python  06_go_enrichment/02_go_enrichment.py
Rscript 06_go_enrichment/03_plot_go_panels.R
```

Step 02 also needs `go-basic.obo` in the results directory:

```bash
curl -o results/06_go_enrichment/go-basic.obo \
     https://release.geneontology.org/2024-01-17/ontology/go-basic.obo
```

## Study sets

Step 00 defines six sets of *O. impluviata* genes from the orthogroups of
section 04:

| Set | Definition |
|---|---|
| `species_specific` | orthogroups found only in *O. impluviata* |
| `cerambycidae_specific` | orthogroups where every species present is a cerambycid |
| `expanded_vs_Aglabripennis` | more genes than in *A. glabripennis* |
| `expanded_vs_Cerambycidae` | more genes than the maximum of the other Cerambycidae |
| `expanded_vs_Chrysomelidae` | more genes than the maximum Chrysomelidae |
| `expanded_vs_Curculionoidea` | more genes than the maximum Curculionoidea |

The reference is the **maximum** of the comparison group rather than its mean,
so an expansion set holds only orthogroups larger in *O. impluviata* than in
every member of the group, and *O. impluviata* is excluded from the "other
Cerambycidae" reference. The thresholds are permissive settings
(`MIN_FOLD_CHANGE=1.0`, `MIN_GENE_DIFF=1`).

## Notes

- The background is every *O. impluviata* gene carrying at least one GO term
  after identifier normalisation, from the eggNOG-mapper annotations. Each
  study set is tested against it with Fisher's exact test, counts propagated up
  the GO hierarchy, and Benjamini–Hochberg correction; terms with FDR < 0.05
  are significant. The three ontologies and the six sets are tested separately.
- Figure 3 labels terms by GO ID because the names do not fit at the panel
  width; step 03 writes the matching `list_GO_terms_*.tsv` alongside each panel.
