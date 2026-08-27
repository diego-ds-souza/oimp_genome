#!/usr/bin/env python3
"""
make_manifest.py

Build the manifest that drives run_collapse_batch.sh and verify_collapse.sh, by
walking an `ncbi datasets download genome` directory tree.

A dehydrated-then-rehydrated datasets download looks like:

    <root>/<anything>/ncbi_dataset/data/<ACCESSION>/protein.faa
    <root>/<anything>/ncbi_dataset/data/<ACCESSION>/genomic.gff
    <root>/<anything>/ncbi_dataset/data/assembly_data_report.jsonl

The species name is taken from assembly_data_report.jsonl when it is present and
from the containing directory name otherwise. Note that for a RefSeq (GCF_)
download the report often names the *paired GenBank* accession and carries no
annotation statistics, so the accession is taken from the data subdirectory --
which is the one the FASTA and GFF actually belong to.

The fifth column, expected_genes, is the number of protein-coding genes in the
GFF. That is the number the collapse must reproduce exactly, so it is what
verify_collapse.sh checks against. It is counted as:

    `gene` features carrying gene_biotype=protein_coding

falling back, for annotations that do not use gene_biotype, to the number of
distinct gene keys among CDS features. The fallback is reported in a comment
line so it is never silently substituted.

Output columns (tab-separated):

    accession  output_tag  protein_fasta  genomic_gff  expected_genes

Usage:
    python3 make_manifest.py NCBI_data/ > manifest.tsv
    python3 make_manifest.py NCBI_data/ --absolute > manifest.tsv
    python3 make_manifest.py NCBI_data/ --full-name > manifest.tsv

Options:
    --absolute    write absolute paths (default: paths as given, so a manifest
                  built with a relative root stays relocatable)
    --full-name   keep the complete organism name, including any subspecies
                  epithet, rather than trimming to the binomial
"""
import argparse
import glob
import gzip
import json
import os
import re
import sys


def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def first_existing(dirpath, *names):
    """Return the first of `names` present in dirpath, plain or gzipped."""
    for n in names:
        for cand in (os.path.join(dirpath, n), os.path.join(dirpath, n + '.gz')):
            if os.path.isfile(cand):
                return cand
    return None


def count_protein_coding_genes(gff):
    """(count, method). Prefers gene_biotype; falls back to CDS gene keys."""
    n_gene = 0
    cds_keys = set()
    with opener(gff) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            ftype, attrs = f[2], f[8]
            if ftype == 'gene' and 'gene_biotype=protein_coding' in attrs:
                n_gene += 1
            elif ftype == 'CDS':
                m = re.search(r'GeneID:(\d+)', attrs)
                if m:
                    cds_keys.add('GeneID:' + m.group(1))
                    continue
                m = re.search(r'(?:^|;)locus_tag=([^;]+)', attrs)
                if m:
                    cds_keys.add('locus_tag:' + m.group(1))
                    continue
                m = re.search(r'(?:^|;)gene=([^;]+)', attrs)
                if m:
                    cds_keys.add('gene:' + m.group(1))
    if n_gene:
        return n_gene, 'gene_biotype=protein_coding'
    return len(cds_keys), 'distinct gene keys among CDS (no gene_biotype in GFF)'


def organism_name(datadir, accession):
    """Species name from assembly_data_report.jsonl, if it is there.

    The report holds one JSON object per assembly in the package, and a package
    can contain more than one -- a RefSeq download often ships the paired
    GenBank assembly alongside it. Prefer the record whose accession matches the
    directory we are actually reading, and fall back to the first record only if
    none matches.
    """
    report = os.path.join(os.path.dirname(datadir), 'assembly_data_report.jsonl')
    if not os.path.isfile(report):
        return None
    first = None
    try:
        with open(report) as fh:
            for line in fh:
                try:
                    d = json.loads(line)
                except ValueError:
                    continue
                name = d.get('organism', {}).get('organismName')
                if not name:
                    continue
                if d.get('accession') == accession:
                    return name
                if first is None:
                    first = name
    except OSError:
        pass
    return first


def tag_from_dirname(dirname):
    """`GCF_000390285.2_Anoplophora_glabripennis` -> `Anoplophora_glabripennis`."""
    s = re.sub(r'^GC[AF]_\d+\.\d+[_\s]*', '', dirname)
    return re.sub(r'[^A-Za-z0-9]+', '_', s).strip('_')


def make_tag(name, full_name):
    parts = re.sub(r'[^A-Za-z0-9]+', ' ', name).split()
    if not full_name:
        parts = parts[:2]
    return '_'.join(parts)


def main():
    ap = argparse.ArgumentParser(
        description='Build a collapse manifest from an ncbi-datasets download.')
    ap.add_argument('root', help='directory containing the per-assembly folders')
    ap.add_argument('--absolute', action='store_true',
                    help='write absolute paths instead of paths as given')
    ap.add_argument('--full-name', action='store_true',
                    help='keep subspecies epithets in the output tag')
    args = ap.parse_args()

    pattern = os.path.join(args.root, '*', 'ncbi_dataset', 'data', '*')
    datadirs = sorted(d for d in glob.glob(pattern) if os.path.isdir(d))
    if not datadirs:
        # also accept a root that *is* an ncbi_dataset/data directory
        datadirs = sorted(d for d in glob.glob(os.path.join(args.root, '*'))
                          if os.path.isdir(d))
    if not datadirs:
        sys.exit(f'make_manifest.py: no assembly directories found under {args.root}')

    rows = []
    notes = []
    skipped = []

    for datadir in datadirs:
        accession = os.path.basename(datadir)
        faa = first_existing(datadir, 'protein.faa')
        gff = first_existing(datadir, 'genomic.gff')
        if not faa or not gff:
            skipped.append((accession, 'no protein.faa' if not faa else 'no genomic.gff'))
            continue

        # the assembly folder is two levels above ncbi_dataset/data/<ACC>
        assembly_dir = os.path.basename(
            os.path.dirname(os.path.dirname(os.path.dirname(datadir))))
        name = organism_name(datadir, accession)
        tag = make_tag(name, args.full_name) if name else tag_from_dirname(assembly_dir)

        n, method = count_protein_coding_genes(gff)
        if method.startswith('distinct'):
            notes.append(f'# NOTE {accession}: expected_genes from {method}')
        if n == 0:
            skipped.append((accession, 'no protein-coding genes found in GFF'))
            continue

        if args.absolute:
            faa, gff = os.path.abspath(faa), os.path.abspath(gff)
        rows.append([accession, tag, faa, gff, str(n)])

    # disambiguate any repeated tag by appending its accession
    seen = {}
    for r in rows:
        seen.setdefault(r[1], []).append(r)
    for tag, group in seen.items():
        if len(group) > 1:
            for r in group:
                r[1] = f'{tag}_{r[0].replace(".", "_")}'
            notes.append(f'# NOTE tag "{tag}" was not unique; accessions appended')

    out = sys.stdout
    out.write('# manifest.tsv  --  input for run_collapse_batch.sh '
              'and verify_collapse.sh\n')
    out.write('# accession\toutput_tag\tprotein_fasta\tgenomic_gff\texpected_genes\n')
    out.write('# expected_genes = protein-coding genes in that GFF; the collapse '
              'must reproduce it exactly\n')
    out.write(f'# {len(rows)} assemblies from {args.root}\n')
    out.write('# fields are TAB-separated and paths may contain spaces: '
              'do not align the columns\n')
    for acc, why in skipped:
        out.write(f'# SKIPPED {acc}: {why}\n')
    for n in notes:
        out.write(n + '\n')
    for a, t, f, g, e in rows:
        out.write(f'{a}\t{t}\t{f}\t{g}\t{e}\n')

    print(f'make_manifest.py: {len(rows)} assemblies written', file=sys.stderr)
    for acc, why in skipped:
        print(f'  SKIPPED {acc}: {why}', file=sys.stderr)
    if skipped:
        # Not an error: a datasets package routinely contains a data directory
        # with only a genome FASTA -- typically the GenBank assembly paired with
        # a RefSeq download. The skips are recorded as comments in the manifest
        # so they are visible there too. Check the list against what you expect.
        print(f'  ({len(skipped)} directory/directories had no proteome; '
              f'recorded as # SKIPPED comments in the manifest)', file=sys.stderr)
    if not rows:
        sys.exit('make_manifest.py: no usable assemblies found')


if __name__ == '__main__':
    main()
