# 02 — Genome assembly and chromosome-scale inference

Methods section 3: contig assembly, haplotig purging, assembly QC, scaffold
naming, telomere identification and synteny (Fig. 1).

## Scripts

| Script | Software | Purpose |
|---|---|---|
| `01_assemble.sh` | hifiasm, any2fasta | Primary contig assembly; graph → FASTA |
| `02_purge_dups.sh` | minimap2, purge_dups | Remove haplotypic duplications |
| `03_assembly_qc.sh` | QUAST, BUSCO | Assembly statistics and completeness |
| `04_sort_and_rename.sh` | seqkit | Sort contigs by size, rename `scaffold_N` |
| `05_fix_contamination.sh` | seqkit, samtools, gffread | Apply the NCBI contamination screen, renumber |
| `06_telomere_repeats.sh` | TIDK | Find and localize telomeric repeats |
| `07_rank_telomeric_repeats.py` | — | Rank the Coleoptera motifs by abundance |
| `08_telomere_summary.py` | — | Call telomeres at scaffold ends |
| `09_pseudokaryotype.py` | matplotlib | Pseudo-karyotype (Fig. 1a) |
| `10_synteny.sh` | minimap2, Circos | Align to *M. alternatus*, bundle links (Fig. 1b) |

## Running

```bash
conda env create -f 02_assembly/environment.yml
conda activate oimp_02_assembly

bash    02_assembly/01_assemble.sh
bash    02_assembly/02_purge_dups.sh
bash    02_assembly/03_assembly_qc.sh
bash    02_assembly/04_sort_and_rename.sh
bash    02_assembly/05_fix_contamination.sh
bash    02_assembly/06_telomere_repeats.sh
python  02_assembly/07_rank_telomeric_repeats.py
python  02_assembly/08_telomere_summary.py
python  02_assembly/09_pseudokaryotype.py
bash    02_assembly/10_synteny.sh
```

## Notes

- Scaffold naming takes two passes. Step 04 sorts the purged contigs by length
  and names them `scaffold_N`; step 05 then applies the NCBI contamination
  screen — removing one scaffold and splitting another at an adaptor match —
  and sorts and renames again, because those edits change the length order.
  Both are needed to reproduce the published scaffold numbers, so run them
  before annotation. Set `GFF_IN` and run step 05 again after annotation to
  bring the GFF3 into the same coordinates.
- Step 05 reproduces edits already applied to GCA_060235215.1; it is here for
  the record.
- Step 08 calls a telomere present when a window within 50 kb of a scaffold end
  carries at least 30 copies of the motif.
- Step 10 aligns to *Monochamus alternatus* (GCA_037114965.1) against the
  unmasked assembly, keeps hits > 100 bp with MAPQ ≥ 30, and bundles links
  within a 1 Mb window with `bundlelinks`.
- Both panels of Fig. 1 were finished in Illustrator; these scripts produce the
  underlying plots.
