# data/

Input files go here. Nothing in this directory is tracked by git.

Download from BioProject
[PRJNA1256903](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1256903):

| File | Accession | Used by |
|---|---|---|
| `hifi_reads.fastq.gz` | SRR38808688 | 01, 02 |
| `rnaseq/raw/<sample>/<sample>_R[12].fastq.gz` | SRR38815773–SRR38815781 | 01, 03, 07 |
| `Oimp_genome.fasta` | GCA_060235215.1 | 02, 03, 08 |
| `Oimp_annotation.gff3` | GCA_060235215.1 | 03, 07, 08 |
| `Oimp_protein.faa` | GCA_060235215.1 | 04 |
| `proteomes/*.faa` | see Table S2 | 04 |

The nine RNA-seq samples are `AF1An`, `AF1Th`, `AF2An`, `AF2Th`, `AM1An`,
`AM2An`, `AM2Th`, `L1He` and `L1Mg`. `01_sequencing/03_rnaseq_fastp.sh` reads
them from `rnaseq/raw/` and writes the trimmed pairs to `rnaseq/trimmed/`,
where sections 03 and 07 look for them.

Rename files as above, or override the path when running a script.
