# 04 — Orthology and phylogeny

Methods section 6, first part: assembling the comparative protein sets,
inferring orthogroups, building the species tree (Fig. 2) and classifying
orthogroups.

## Scripts

| Script | Software | Purpose |
|---|---|---|
| `00_download_proteomes.sh` | NCBI datasets | Fetch the 24 comparison proteomes and their GFF3s |
| `01_busco_proteomes.sh` | compleasm | Completeness of each proteome as downloaded |
| `02_collapse_isoforms.sh` | — | Reduce each proteome to one protein per gene |
| `03_run_orthofinder.sh` | OrthoFinder, DIAMOND | Orthogroup inference |
| `04_align_and_trim_sco.sh` | MAFFT, trimAl | Align and trim the single-copy orthologs |
| `05_concat_and_partitions.sh` | — | Supermatrix and partition file |
| `06_iqtree.sh` | IQ-TREE | Partitioned ML tree with support values |
| `07_root_tree.R` | ape | Root on *T. castaneum*, ladderize, write the taxonomy table |
| `08_orthology_categories.sh` | — | Classify orthogroups into the nine categories |
| `09_verify_counts.sh` | — | Cross-check the category counts against OrthoFinder |
| `10_orthology_histogram.R` | ggplot2 | Orthology distribution figure |

`isoform_collapse/` is the toolkit called by step 02 and has its own README.

## Running

```bash
conda env create -f 04_orthology/environment.yml
conda activate oimp_04_orthology

bash    04_orthology/00_download_proteomes.sh
bash    04_orthology/01_busco_proteomes.sh
bash    04_orthology/02_collapse_isoforms.sh
bash    04_orthology/03_run_orthofinder.sh
bash    04_orthology/04_align_and_trim_sco.sh
bash    04_orthology/05_concat_and_partitions.sh
bash    04_orthology/06_iqtree.sh
Rscript 04_orthology/07_root_tree.R
bash    04_orthology/08_orthology_categories.sh
bash    04_orthology/09_verify_counts.sh
Rscript 04_orthology/10_orthology_histogram.R
```

Steps 03 onwards locate the OrthoFinder output themselves, so a re-run needs no
path edited. OrthoFinder in MSA mode is the expensive step.

## Notes

- The panel is *O. impluviata* plus 24 annotated beetle genomes covering
  Cerambycidae, Chrysomelidae and Curculionoidea, with *Tribolium castaneum* as
  the outgroup. Accessions are in `accessions.txt`.
- Comparative methods assume one protein per gene. NCBI protein FASTA headers
  carry no gene identifier, so step 02 maps each protein to its parent gene
  through the genomic GFF3 rather than through header text, and
  `isoform_collapse/verify_collapse.sh` checks the result.
- Completeness was measured with compleasm against the insecta_odb10 BUSCO
  marker set, on the proteomes **as downloaded**, before isoform reduction.
- IQ-TREE ran with `-m MFP+MERGE`, 1,000 ultrafast bootstrap replicates and
  1,000 SH-aLRT tests. Step 07 roots the tree and writes `species_taxonomy.tsv`,
  the family assignments used by step 08.
- Step 08 sorts orthogroups into single-copy universal, multicopy universal,
  species-specific, Cerambycidae-, Chrysomelidae-, Curculionoidea- and
  Phytophaga-specific, other shared orthologs, and unassigned genes. Step 09
  recomputes the totals from `Orthogroups.GeneCount.tsv` as an independent
  check.
- The gene family **expansion** sets are defined in
  `06_go_enrichment/00_build_gene_sets.sh`, where they are used.
