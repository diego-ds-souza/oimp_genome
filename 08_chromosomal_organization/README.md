# 08 — Chromosomal organization of host-interface genes

Methods section 8: placing the host-interface gene families on chromosome-scale
coordinates, finding tandem clusters and superloci, and testing whether
differentially expressed genes are concentrated in those clusters (Fig. 7).

## Scripts

| Script | Software | Purpose |
|---|---|---|
| `01_chromosome_ideograms.R` | ggplot2 | Chromosome coordinates, density, clusters, superloci, ideograms |
| `02_build_cluster_tables.py` | openpyxl | Cluster and superlocus membership tables |
| `03_cluster_de_enrichment.R` | — | Fisher and Mantel–Haenszel tests of DE enrichment in clusters |

## Running

```bash
conda env create -f 08_chromosomal_organization/environment.yml
conda activate oimp_08_chromosomal_organization

Rscript 08_chromosomal_organization/01_chromosome_ideograms.R
python  08_chromosomal_organization/02_build_cluster_tables.py
Rscript 08_chromosomal_organization/03_cluster_de_enrichment.R
```

`scaffolds_to_chromosome.csv` maps each scaffold to a chromosome, records the
telomere calls from section 02 and whether the synteny supported the
assignment; step 01 reads it.

## Definitions

- A **tandem cluster** is two or more genes of the same family on the same
  chromosome, chained while consecutive genes are within 100 kb.
- A **superlocus** is a run of clusters merged while each starts within 50 kb of
  the running end, counted only when the run spans more than one family.
- Gene density is computed in 500 kb sliding windows.

All four thresholds are settings in `01_chromosome_ideograms.R`.

## Notes

- Gene positions in the GFF3 are scaffold-relative. Step 01 converts them by
  concatenating the scaffolds assigned to each chromosome in order. Scaffold
  orientation is set from the telomere positions: where the telomeric arrays
  imply the opposite polarity the scaffold is reversed before its coordinates
  are made cumulative, so every chromosome runs in a consistent direction.
  Ideograms mark those junctions with a double diagonal.
- Step 01 joins one gene list per family onto the gene positions, so a gene
  belonging to two families is clustered twice, once under each. Step 02
  collapses each gene to a single primary family first (OR over GR; UGT and
  SerPro over CAZyme), re-derives the clusters and superloci with the same
  rules, and checks its own counts against the published table, printing any
  cluster that disagrees.
- Every clustered gene is by construction a member of one of the 16 families,
  and those families are differentially expressed more often than the average
  annotated gene, so a genome-wide 2×2 table compares groups differing in
  family membership as well as in clustering. Step 03 runs that genome-wide
  test, repeats it restricted to the 16 families, runs it separately for
  stage-biased and for tissue- or sex-biased genes, and reports a
  Mantel–Haenszel estimate stratified by family. The odds ratio reported is the
  sample odds ratio, (a·d)/(b·c), not the conditional maximum-likelihood
  estimate `fisher.test()` returns.
