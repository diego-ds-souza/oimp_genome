# 03 — Gene prediction and functional annotation

Methods section 5: repeat identification and masking, the Funannotate
gene-prediction and annotation workflow, curation of the final gene set, and
its completeness.

## Scripts

| Script | Software | Purpose |
|---|---|---|
| `01_repeat_library.sh` | RepeatModeler | De novo repeat families; custom library |
| `02_repeatmasker.sh` | RepeatMasker | Annotate and soft-mask repeats; divergence summary |
| `03_divsum_table.py` | — | Extract the class-by-divergence matrix |
| `04_repeat_landscape.R` | ggplot2 | Repeat landscape (Fig. S1b) |
| `05_funannotate_train.sh` | Funannotate | Genome-guided assembly (HISAT2, StringTie), GMAP/BLAT, PASA |
| `06_funannotate_predict.sh` | Funannotate | Ab initio predictors, Exonerate, EVidenceModeler |
| `07_funannotate_update.sh` | Funannotate | PASA refinement, UTRs, kallisto quantification |
| `08_interproscan.sh` | InterProScan | Domain and GO assignment |
| `09_funannotate_annotate.sh` | Funannotate | Pfam, eggNOG/COG, UniProt, MEROPS, dbCAN, GO |
| `10_gene_evidence.py` | — | Score each gene model by its supporting evidence |
| `11_filter_abinitio_only.py` | — | Drop models resting on ab initio prediction alone |
| `12_extract_curated_annotation.sh` | AGAT, seqkit | Subset the GFF3 and proteins to the curated set |
| `13_annotation_busco.sh` | BUSCO | Completeness of the protein set |

`helper.functions.R`, used by `04_repeat_landscape.R`, was written by
Terrence Sylvester and is included unchanged.

## Running

```bash
conda env create -f 03_annotation/environment.yml
conda activate oimp_03_annotation

bash    03_annotation/01_repeat_library.sh
bash    03_annotation/02_repeatmasker.sh
python  03_annotation/03_divsum_table.py
Rscript 03_annotation/04_repeat_landscape.R
bash    03_annotation/05_funannotate_train.sh
bash    03_annotation/06_funannotate_predict.sh
bash    03_annotation/07_funannotate_update.sh
bash    03_annotation/08_interproscan.sh
bash    03_annotation/09_funannotate_annotate.sh
python  03_annotation/10_gene_evidence.py
python  03_annotation/11_filter_abinitio_only.py
bash    03_annotation/12_extract_curated_annotation.sh
bash    03_annotation/13_annotation_busco.sh
```

Steps 05–09 share one output directory and must run in order. Funannotate needs
its databases installed first (`funannotate setup -d $FUNANNOTATE_DB -b
insecta`) and GeneMark-ES needs a licence key.

## Notes

- Step 05 reads the trimmed RNA-seq pairs from
  [`01_sequencing/03_rnaseq_fastp.sh`](../01_sequencing/README.md); the
  sample-specific Trinity assemblies from `01_sequencing/04_trinity_denovo.sh`
  were supplied to it as additional transcript evidence.
- The de novo repeat families were combined with curated Coleoptera repeats
  from Dfam 3.8 and Repbase, and RepeatMasker was run against that library in
  its most sensitive mode (`-s`). Masking for gene prediction was an earlier,
  separate run on the pre-renumbering assembly; the repeat content reported in
  the paper comes from the run in step 02, on the final assembly.
- Step 10 scores each model on three kinds of support: **protein** — an
  InterPro or Pfam cross-reference, an eggNOG, COG or MEROPS note, or a GO
  term; **transcript** — 5′ or 3′ UTR features, which Funannotate adds only
  where PASA matched an assembled transcript; and **ab initio**, true of every
  model. Step 11 removes the models with ab initio support alone, with no
  rescue rule.
- The curated GFF3 was standardised with AGAT and given the final scaffold
  coordinates by `02_assembly/05_fix_contamination.sh`.
