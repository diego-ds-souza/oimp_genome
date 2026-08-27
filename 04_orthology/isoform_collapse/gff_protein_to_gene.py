#!/usr/bin/env python3
"""
gff_protein_to_gene.py

Emit a two-column map, protein accession -> gene key, from an NCBI genomic.gff.

NCBI `protein.faa` headers carry no gene identifier, so the link between a
protein and its gene has to be read from the annotation. Every CDS line in
genomic.gff carries `protein_id=` alongside a gene identifier:

    CDS ... Dbxref=GeneID:108903154,GenBank:XP_018560752.1;gene=LOC108903154;
            protein_id=XP_018560752.1

Gene key priority, first match wins:
    1. Dbxref=GeneID:<n>   stable across releases; present in RefSeq
    2. locus_tag=<tag>     stable; present in most GenBank submissions
    3. gene=<symbol>       a symbol, not guaranteed unique, so used last
    4. the Parent chain    CDS -> mRNA -> gene feature

Protein accessions are written WITHOUT their version suffix: `XP_018560752.1`
becomes `XP_018560752`. collapse_isoforms_from_gff.py strips the version the
same way when it writes `>{tag}|{accession}`, so the join key matches the
identifiers in an OrthoFinder run built from those files -- but it is not a
literal NCBI accession, and joining against anything that keeps the version
suffix needs that suffix stripped first.

Usage:
    python3 gff_protein_to_gene.py --gff ncbi/GCF_000390285.2/genomic.gff \\
                                   --out maps/GCF_000390285.2.prot2gene.tsv

Requires Python 3.6 or newer. Standard library only. Gzipped input is accepted.
"""
import argparse
import gzip
import re


def opener(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def parse_gff(path):
    prot2gene, mrna2gene, gene_of_id, pending = {}, {}, {}, []

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

            if ftype in ('gene', 'pseudogene'):
                gid, k = attr('ID'), gene_key()
                if gid and k:
                    gene_of_id[gid] = k
            elif ftype in ('mRNA', 'transcript'):
                tid, k = attr('ID'), gene_key()
                if tid:
                    mrna2gene[tid] = k if k else ('PARENT', attr('Parent'))
            elif ftype == 'CDS':
                pid = attr('protein_id')
                if not pid:
                    continue
                pid = re.sub(r'\.\d+$', '', pid)
                if pid in prot2gene:
                    continue                      # multi-exon CDS, already seen
                k = gene_key()
                if k:
                    prot2gene[pid] = k
                else:
                    pending.append((pid, attr('Parent')))

    for pid, parent in pending:
        if pid in prot2gene or not parent:
            continue
        k = mrna2gene.get(parent)
        if isinstance(k, tuple) and k[0] == 'PARENT':
            k = gene_of_id.get(k[1])
        if k:
            prot2gene[pid] = k
    return prot2gene


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gff', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    m = parse_gff(args.gff)
    with open(args.out, 'w') as fh:
        fh.write('protein_accession\tgene_key\n')
        for p in sorted(m):
            fh.write(f'{p}\t{m[p]}\n')

    genes = len(set(m.values()))
    print(f'{args.gff}')
    print(f'  proteins mapped : {len(m):,}')
    print(f'  distinct genes  : {genes:,}')
    print(f'  isoforms        : {len(m) - genes:,}')
    print(f'  written         : {args.out}')


if __name__ == '__main__':
    main()
