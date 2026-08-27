# 01 — Sequencing and read processing

Methods sections 2 and 4: read statistics and genome-size estimation for the
PacBio HiFi reads, and trimming and de novo assembly of the RNA-seq libraries.

## Scripts

| Script | Software | Purpose |
|---|---|---|
| `01_read_stats.sh` | seqkit | Read count, total bases, N50, quality, GC |
| `02_kmer_count.sh` | Jellyfish | Canonical 21-mer counts and histogram |
| `03_rnaseq_fastp.sh` | fastp | Adapter and quality trimming of the RNA-seq libraries |
| `04_trinity_denovo.sh` | Trinity | Sample-specific de novo transcriptome assemblies |

## Running

```bash
conda env create -f 01_sequencing/environment.yml
conda activate oimp_01_sequencing

bash 01_sequencing/01_read_stats.sh
bash 01_sequencing/02_kmer_count.sh
bash 01_sequencing/03_rnaseq_fastp.sh
bash 01_sequencing/04_trinity_denovo.sh
```

Steps 01–02 act on the HiFi reads, steps 03–04 on the RNA-seq libraries.

## Notes

- Genome size was estimated by uploading the histogram from step 02 to the
  [GenomeScope 2.0 server](http://genomescope.org/genomescope2.0/) with k = 21
  and the other settings left at their defaults. `genomescope2` run locally
  gives the same output.
- Jellyfish cannot read gzipped FASTQ; decompress first, or pass
  `READS=<(zcat data/hifi_reads.fastq.gz)`.
- Step 03 writes trimmed pairs, merged pairs, unpaired reads and per-library
  reports to `data/rnaseq/trimmed/<sample>/`. The **unmerged pairs** are the
  input to sections 03 and 07.
