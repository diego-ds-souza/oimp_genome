#!/usr/bin/env python3
"""
collapse_isoforms_from_gff.py

Reduce an NCBI proteome to one protein per gene, using the annotation GFF
rather than the text of the FASTA header.

Why the GFF. NCBI `protein.faa` headers carry no gene identifier:

    >XP_018560752.1 uncharacterized protein LOC108903154 [Anoplophora glabripennis]
    >NP_001034280.2 dorsocross [Tribolium castaneum]
    >KAJ8931856.1 hypothetical protein NQ318_009800, partial [Aromia moschata]

There is no `gene=` and no `locus_tag=` field, so any rule that looks for one
falls through to the protein accession, which is unique to each isoform and
therefore collapses nothing. `genomic.gff` does carry the link: every CDS line
has `protein_id=` together with `Dbxref=GeneID:`, `gene=` or `locus_tag=`.

Gene key priority (first match wins):
    1. Dbxref=GeneID:<n>     stable across releases, present in RefSeq
    2. locus_tag=<tag>       stable, present in most GenBank submissions
    3. gene=<symbol>         a symbol, not always unique, so used only as a fallback
    4. Parent gene feature ID resolved through the mRNA

For each gene the longest protein is kept. Ties on length are broken by keeping
the lexicographically largest accession, so the result is deterministic and does
not depend on the order of the input file.

Output headers are `>{tag}|{accession}` with the version suffix stripped
(`XP_066262389.1` is written as `XP_066262389`). That is the same key
gff_protein_to_gene.py emits, so the two join directly; it is not, however, a
literal NCBI accession any more.

Verification. When --expected-genes is supplied the script compares the number
of genes recovered against it and sets its exit status:

    exact match          OK        exit 0
    within +/-5%         WARN      exit 0
    outside +/-5%        FAIL      exit 1

With an NCBI GFF the answer should be exact, because the expected value is the
protein-coding gene count of the same annotation the proteins came from. Any
shortfall is reported as `proteins with no GFF gene`.

Usage:
    python3 collapse_isoforms_from_gff.py \
        --fasta  NCBI_data/GCF_000390285.2_.../protein.faa \
        --gff    NCBI_data/GCF_000390285.2_.../genomic.gff \
        --out    clean_proteomes/Anoplophora_glabripennis.fa \
        --tag    Anoplophora_glabripennis \
        --expected-genes 14828

Requires Python 3.6 or newer. Standard library only. Gzipped input is accepted.
"""
import argparse
import gzip
import re
import sys


def opener(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def parse_gff(path):
    """protein accession -> gene key, from the CDS features of an NCBI GFF."""
    prot2gene = {}
    mrna2gene = {}          # transcript ID -> gene key, for the Parent fallback
    gene_of_id = {}         # gene feature ID -> gene key
    pending = []            # CDS whose gene key must come from the Parent chain

    with opener(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            ftype, attrs = f[2], f[8]

            def attr(name):
                m = re.search(rf'(?:^|;){name}=([^;]+)', attrs)
                return m.group(1) if m else None

            def gene_key():
                d = attr('Dbxref') or attr('dbxref') or ''
                m = re.search(r'GeneID:(\d+)', d)
                if m:
                    return f'GeneID:{m.group(1)}'
                lt = attr('locus_tag')
                if lt:
                    return f'locus_tag:{lt}'
                g = attr('gene')
                if g:
                    return f'gene:{g}'
                return None

            if ftype == 'gene' or ftype == 'pseudogene':
                gid = attr('ID')
                k = gene_key()
                if gid and k:
                    gene_of_id[gid] = k

            elif ftype in ('mRNA', 'transcript'):
                tid = attr('ID')
                k = gene_key()
                if tid:
                    if k:
                        mrna2gene[tid] = k
                    else:
                        p = attr('Parent')
                        if p:
                            mrna2gene[tid] = ('PARENT', p)

            elif ftype == 'CDS':
                pid = attr('protein_id')
                if not pid:
                    continue
                pid = re.sub(r'\.\d+$', '', pid)      # drop the version
                if pid in prot2gene:
                    continue                          # multi-exon CDS, already seen
                k = gene_key()
                if k:
                    prot2gene[pid] = k
                else:
                    pending.append((pid, attr('Parent')))

    # resolve CDS -> mRNA -> gene
    for pid, parent in pending:
        if pid in prot2gene or not parent:
            continue
        k = mrna2gene.get(parent)
        if isinstance(k, tuple) and k[0] == 'PARENT':
            k = gene_of_id.get(k[1])
        if k:
            prot2gene[pid] = k
        else:
            prot2gene[pid] = f'unmapped:{pid}'
    return prot2gene


def read_fasta(path):
    hdr, seq = None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if hdr is not None:
                    yield hdr, ''.join(seq)
                hdr, seq = line.strip(), []
            else:
                seq.append(line.strip())
    if hdr is not None:
        yield hdr, ''.join(seq)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fasta', required=True)
    ap.add_argument('--gff', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--tag', required=True,
                    help='species tag written into the output headers')
    ap.add_argument('--expected-genes', type=int,
                    help='protein-coding gene count for this annotation; with an '
                         'NCBI GFF the collapse should reproduce it exactly')
    ap.add_argument('--min-length', type=int, default=0, metavar='N',
                    help='drop retained proteins shorter than N residues. Applied '
                         'after the collapse, so it never changes which isoform '
                         'is chosen. Counted separately from the gene check.')
    args = ap.parse_args()

    prot2gene = parse_gff(args.gff)
    print(f'GFF          : {args.gff}')
    print(f'  protein accessions mapped : {len(prot2gene):,}')

    best, unmapped = {}, 0
    n_in = 0
    for hdr, seq in read_fasta(args.fasta):
        n_in += 1
        acc = re.sub(r'\.\d+$', '', hdr.lstrip('>').split()[0])
        gene = prot2gene.get(acc)
        if gene is None or gene.startswith('unmapped:'):
            unmapped += 1
            gene = f'unmapped:{acc}'          # keep it, as its own gene
        cur = best.get(gene)
        if cur is None or (len(seq), acc) > (len(cur[1]), cur[0]):
            best[gene] = (acc, seq)

    n_genes = len(best)
    dropped = 0
    with open(args.out, 'w') as fo:
        for gene in sorted(best):
            acc, seq = best[gene]
            if len(seq) < args.min_length:
                dropped += 1
                continue
            fo.write(f'>{args.tag}|{acc}\n')
            for i in range(0, len(seq), 60):
                fo.write(seq[i:i + 60] + '\n')

    print(f'FASTA        : {args.fasta}')
    print(f'  sequences in              : {n_in:,}')
    print(f'  genes out                 : {n_genes:,}')
    print(f'  isoforms collapsed        : {n_in - n_genes:,}')
    print(f'  proteins with no GFF gene : {unmapped:,}')
    if args.min_length:
        label = f'dropped under {args.min_length} aa'
        print(f'  {label:<26}: {dropped:,}')
        print(f'  {"sequences written":<26}: {n_genes - dropped:,}')
    print(f'  written                   : {args.out}')

    status = 0
    if args.expected_genes:
        exp = args.expected_genes
        d = n_genes - exp
        pct = 100.0 * n_genes / exp
        if d == 0:
            flag = 'OK'
        elif 95.0 <= pct <= 105.0:
            flag = '*** WARN ***'
        else:
            flag = '*** FAIL ***'
            status = 1
        print(f'  protein-coding genes      : {exp:,}  '
              f'({pct:.1f}% of expected, {d:+,})  {flag}')
        if d and unmapped:
            print(f'  -> {unmapped:,} proteins had no gene in this GFF, which '
                  f'accounts for some or all of the difference.')
        elif d:
            print('  -> counts differ with no unmapped proteins: check that the '
                  'GFF and the FASTA are from the same assembly and release.')
    elif n_in == n_genes:
        # Only meaningful without an expected count: with one, an exact match
        # already proves the collapse was correct even when it removed nothing.
        print('  NOTE: nothing was collapsed. That is correct for an annotation '
              'with one protein per gene; otherwise confirm the GFF matches this '
              'proteome and that its CDS lines carry protein_id=.', file=sys.stderr)

    return status


if __name__ == '__main__':
    sys.exit(main())
