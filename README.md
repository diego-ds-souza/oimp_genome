# *Oncideres impluviata* genome — analysis code

Scripts used in **"Chromosome-scale genome of a twig-girdler longhorn beetle
reveals clustered genomic organization of herbivory-linked functions"**
(Souza et al., *BMC Genomics*).

Sequence data and large intermediate files are not in this repository; see
[Data](#data).

## Repository layout

Directories follow the order of the Methods.

| Directory | Methods section |
|---|---|
| `01_sequencing/` | 2. DNA extraction and genome sequencing; 4. RNA extraction and sequencing |
| `02_assembly/` | 3. Genome assembly and chromosome-scale inference |
| `03_annotation/` | 5. Gene prediction and annotation |
| `04_orthology/` | 6. Comparative genomics — orthology and phylogeny |
| `05_gene_family_evolution/` | 6. Comparative genomics — gene family evolution |
| `06_go_enrichment/` | 6. Comparative genomics — GO enrichment |
| `07_transcriptomics/` | 7. Transcriptomic profiling and gene family analysis |
| `08_chromosomal_organization/` | 8. Chromosomal organization of host-interface genes |

Each directory holds its scripts in run order and a README listing them.

## Running

Scripts are run **from the repository root**, with no arguments:

```bash
conda env create -f 01_sequencing/environment.yml
conda activate oimp_01_sequencing
bash 01_sequencing/02_kmer_count.sh
```

Inputs are read from `data/` and results written to `results/<section>/`;
neither directory is tracked by git. See [`data/README.md`](data/README.md).

Paths, thread counts and parameters are in a settings block at the top of each
script. Edit them, or override any of them on the command line:

```bash
READS=data/my_reads.fastq THREADS=32 bash 01_sequencing/02_kmer_count.sh
```

The scripts in `04_orthology/isoform_collapse/` are general-purpose tools and
take command-line arguments instead.

Each section has its own `environment.yml`, named `oimp_<section>`. R scripts
were run under R 4.3. The sections are meant to be run in order.

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

## Data

| Resource | Accession |
|---|---|
| BioProject | [PRJNA1256903](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1256903) |
| Genome assembly and annotation | [GCA_060235215.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_060235215.1/) |
| PacBio HiFi reads | SRR38808688 |
| RNA-seq libraries | SRR38815773–SRR38815781 |
| Comparison proteomes | see Table S2 |

Result tables are published with the paper as Supplementary Tables S1–S18 and
Supplementary files 1–7. Intermediate files that exceed GitHub's size limits —
assembly graphs, alignments, the OrthoFinder working directory, per-orthogroup
alignments and the raw CAFE output — are omitted and can be regenerated with
these scripts.

## Citation

> Souza D, Sylvester T, et al. Chromosome-scale genome of a twig-girdler
> longhorn beetle reveals clustered genomic organization of herbivory-linked
> functions. *BMC Genomics*. 2026.

## Contact

Diego Souza — tsepulveda@fieldmuseum.org

Released under the MIT License (see `LICENSE`).
