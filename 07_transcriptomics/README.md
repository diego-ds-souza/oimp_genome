# 07 — Transcriptomic profiling and gene family expression

Methods section 7: quantification of the nine RNA-seq libraries, differential
expression, and expression profiles of the host-interface gene families
(Fig. 6).

## Scripts

| Script | Software | Purpose |
|---|---|---|
| `01_parse_gff3_annotations.sh` | — | tx2gene, gene names and family membership from the GFF3 |
| `02_salmon_index.sh` | Salmon | Decoy-aware index |
| `03_salmon_quant.sh` | Salmon | Quantify every library |
| `04_deseq2_analysis.R` | tximport, DESeq2 | Seven differential expression contrasts |
| `05_expression_figures.R` | ggplot2, pheatmap | Family expression, DE summaries, PCA |
| `06_volcano_target_families.R` | ggplot2 | Volcano plots of the target families (Fig. 6) |

## Running

```bash
conda env create -f 07_transcriptomics/environment.yml
conda activate oimp_07_transcriptomics

bash    07_transcriptomics/01_parse_gff3_annotations.sh
bash    07_transcriptomics/02_salmon_index.sh
bash    07_transcriptomics/03_salmon_quant.sh
Rscript 07_transcriptomics/04_deseq2_analysis.R
Rscript 07_transcriptomics/05_expression_figures.R
Rscript 07_transcriptomics/06_volcano_target_families.R
```

`samples.tsv` lists the nine libraries with their stage, sex and tissue, and
the paths to the trimmed FASTQ pairs written by
[`01_sequencing/03_rnaseq_fastp.sh`](../01_sequencing/README.md).

## Notes

- The Salmon index is decoy-aware: the genome is supplied alongside the
  transcripts so reads from genomic regions resembling a transcript are not
  assigned to it. Quantification uses automatic library-type detection with
  sequence-specific, GC and positional bias correction, and abundances are
  summarised to gene level with tximport.
- Differential expression uses size-factor normalisation, the Wald test and
  Benjamini–Hochberg correction. A gene is called differentially expressed at
  adjusted **p < 0.05 and |log2FC| > 1**; the fold-change threshold is applied
  in steps 05 and 06, since step 04 writes the unfiltered DESeq2 tables.
- Seven contrasts are run and they are not equivalent. The larval libraries
  come from a single individual, so the three larva-versus-adult contrasts have
  no biological replication and are descriptive only; the male-versus-female
  thorax contrast has one replicate on one side. Each output file records which
  kind it is.
- Family membership is read from the mRNA attributes of the curated annotation
  — diagnostic Pfam domains, MEROPS classifications and CAZy assignments — so
  it matches what was deposited rather than a separate curation.
- Family-level expression is the **mean** normalised expression across the
  genes of a family, which keeps families of very different sizes comparable.
