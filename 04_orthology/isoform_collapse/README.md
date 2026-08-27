# isoform_collapse — one protein per gene

Reduces downloaded NCBI proteomes to one protein per gene, which is what
OrthoFinder, CAFE and gene family counting assume they are given. Called by
`04_orthology/02_collapse_isoforms.sh`; the scripts here are general-purpose
and take command-line arguments.

## Why not the FASTA headers

NCBI `protein.faa` headers carry no gene identifier — no `gene=`, no
`locus_tag=`. A script that keys on those fields falls back to the protein
accession, which is unique per isoform, so it collapses nothing and exits
successfully. Grouping by description after stripping ` isoform X` fails the
other way, because thousands of proteins share `hypothetical protein`.
The ` isoform X` convention is RefSeq's; GenBank proteomes carry redundancy
that no header rule can see. OrthoFinder's own `tools/primary_transcript.py`
looks for Ensembl-style `gene:` tags and falls through in the same way.

The rule is therefore: **take the gene key from the annotation, never from the
FASTA header.** Every CDS line in `genomic.gff` links protein to gene. The key
is the first available of `Dbxref=GeneID:`, `locus_tag=`, `gene=`, or the
CDS → mRNA → gene `Parent` chain; the longest protein per gene is kept, ties
broken by accession so the result is deterministic. Output headers are
`>{tag}|{accession}` with the version suffix stripped.

Sequence-identity clustering (CD-HIT, MMseqs2) is not a substitute: it merges
recent duplicates — the tandem families a comparative analysis is there to
count — and misses divergent isoforms of one gene.

## Usage

```bash
python3 make_manifest.py NCBI_data/ > manifest.tsv
bash run_collapse_batch.sh manifest.tsv clean_proteomes/
bash verify_collapse.sh   manifest.tsv clean_proteomes/
orthofinder -f clean_proteomes/ -t 16
```

Or one genome at a time:

```bash
python3 collapse_isoforms_from_gff.py \
    --fasta protein.faa --gff genomic.gff \
    --out clean_proteomes/Anoplophora_glabripennis.fa \
    --tag Anoplophora_glabripennis --expected-genes 14828
```

## The manifest

Tab-separated, one row per genome:

```
accession <TAB> output_tag <TAB> protein_fasta <TAB> genomic_gff <TAB> expected_genes
```

`make_manifest.py` builds it from an `ncbi datasets` download tree, taking the
species name from `assembly_data_report.jsonl` and counting `expected_genes`
from each GFF. Fields must be TAB-separated, because download folder names
contain spaces, and the file must have Unix line endings.

Write `-` in the GFF field for a genome that arrives already collapsed; the
batch driver then leaves the file alone and the verifier checks only that it
holds the number of sequences claimed.

Outputs are named `<output_tag>.fa`, with a single extension, because
OrthoFinder takes the species name from the filename with one extension
stripped.

## Verification

**The test is arithmetic: the collapsed sequence count must equal the number of
protein-coding genes in the annotation the proteins came from** — equal, not
close. Count those genes as:

```bash
awk -F'\t' '$3=="gene" && $9 ~ /gene_biotype=protein_coding/' genomic.gff | wc -l
```

Not the total number of `gene` features, which includes lncRNA, tRNA, rRNA and
pseudogene loci that contribute no protein, and against which a correct run
looks like it lost a fifth of its genes.

`verify_collapse.sh` applies the test per genome and also reports any proteome
present in the clean directory but absent from the manifest, since OrthoFinder
reads the whole directory and such a file would otherwise be analysed
unchecked. `collapse_isoforms_from_gff.py --expected-genes` sets its exit
status the same way: exact `OK`, within ±5% `WARN`, outside `FAIL`.

"Nothing was collapsed" is sometimes correct — some annotations already have
one protein per gene — and an exact match proves it.

## Files

| File | Purpose |
|---|---|
| `make_manifest.py` | download tree → `manifest.tsv`, with gene counts |
| `collapse_isoforms_from_gff.py` | one proteome → one protein per gene |
| `gff_protein_to_gene.py` | GFF → `protein_accession <TAB> gene_key` |
| `run_collapse_batch.sh` | batch driver over a manifest |
| `verify_collapse.sh` | manifest-driven check, raw vs collapsed vs expected |

The Python scripts need Python 3.6+ and the standard library only, and accept
gzipped input. The shell scripts avoid bash-4-only syntax.
