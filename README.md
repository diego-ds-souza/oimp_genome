# *Oncideres impluviata* genome — analysis code

This repository holds the scripts used in **"Chromosome-scale genome of a twig-girdler longhorn beetle
reveals clustered genomic organization of herbivory-linked functions"**
(Souza et al., *BMC Genomics*, under review).

## Repository

Directories follow the order of the Methods in the paper.

| Directory | Methods section | Contents |
|---|---|---|
| `01_sequencing/` | 2. DNA extraction and genome sequencing; 4. RNA extraction and sequencing | Read statistics; *k*-mer counting and genome-size estimation; RNA-seq trimming and de novo transcriptome assembly |
| `02_assembly/` | 3. Genome assembly and chromosome-scale inference | Assembly, haplotig purging, assembly QC, scaffold naming, telomere identification, synteny |
| `03_annotation/` | 5. Gene prediction and annotation | Repeat modelling and masking; structural and functional annotation; annotation completeness |
| `04_orthology/` | 6. Comparative genomics | Proteome retrieval, isoform collapsing, orthology inference, phylogeny, orthogroup classification |
| `05_gene_family_evolution/` | 6. Comparative genomics | Birth–death modelling of gene family size with CAFE; host-interface family summary |
| `06_go_enrichment/` | 6. Comparative genomics | Study gene sets, GO mapping, enrichment of expanded orthogroups |
| `07_transcriptomics/` | 7. Transcriptomic profiling and gene family analysis | Family assignment, quantification, differential expression, expression figures |
| `08_chromosomal_organization/` | 8. Chromosomal organization of host-interface genes | Chromosome coordinates, tandem clusters, superloci, ideograms, enrichment tests |

---

Each directory holds its scripts in run order and a README listing them.

## Running

Scripts are run **from the repository root**, with no arguments:

```bash
conda activate oimp_01_sequencing
bash 01_sequencing/02_kmer_count.sh
```

**Paths are relative to the repository root.** Inputs are read from `data/`
and results are written to `results/<section>/`. 

**Settings live in a block at the top of each script**, so the file records
exactly what was run:

```bash
# ------------------------------- settings ------------------------------------
READS="${READS:-data/hifi_reads.fastq}"        # PacBio HiFi reads, uncompressed
K="${K:-21}"                                   # k-mer length (21 for GenomeScope)
THREADS="${THREADS:-10}"                       # CPU threads
# -----------------------------------------------------------------------------
```

Edit them, or override any of them on the command line:

```bash
READS=data/my_reads.fastq THREADS=32 bash 01_sequencing/02_kmer_count.sh
```

The few scripts that are general-purpose tools rather than a record of one
analysis take command-line arguments instead; each says so in its header and
prints usage when run with `-h`.

**Environments.** Each section has its own `environment.yml`, named
`oimp_<section>`:

```bash
conda env create -f 01_sequencing/environment.yml
conda activate oimp_01_sequencing
```

R scripts were run under R 4.3 with `ape`, `ggplot2`, `ggtree`, `DESeq2`,
`circlize` and `RColorBrewer`. 

The sections are meant to be run in order.

## Software

| Step | Software |
|---|---|
| Genome size | Jellyfish v2.3.1, GenomeScope v2.0 |
| Read processing | fastp v1.0.1, Trinity |
| Assembly | Hifiasm v0.16.1, any2fasta v0.4.2, Purge_dups v1.2.6 |
| Assembly QC | QUAST v5.0.2, BUSCO v5.2.2 (insecta_odb10) |
| Chromosome inference | TIDK v0.2.65, minimap2, Circos v0.23, circlize v0.4.17 |
| Repeats | RepeatModeler v2.0.6, RepeatMasker v4.1.8, rmblast v2.14.1, Dfam 3.8 |
| Annotation | Funannotate v1.8.17, InterProScan, eggNOG-mapper, dbCAN, MEROPS |
| Orthology | OrthoFinder v2.5.5, DIAMOND, compleasm |
| Phylogeny | MAFFT v7.520, trimAl v1.4, AMAS v1.0, IQ-TREE v3.1.3, ape 5.0 |
| Gene family evolution | CAFE v5.1.0 |
| GO enrichment | GOATOOLS v1.3.1, go-basic.obo 2024-01-17 |
| Expression | Salmon v1.10, DESeq2 v1.38 |

---

## Citation

> Souza DS, et al. Chromosome-scale genome of a twig-girdler
> longhorn beetle reveals clustered genomic organization of herbivory-linked
> functions. *BMC Genomics*. 2026 (under review).
