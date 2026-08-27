#!/usr/bin/env python3
"""
Build the CAFE5 gene family count table from the OrthoFinder output.

CAFE5 expects a tab-delimited table whose first two columns are a description
and a family identifier, followed by one integer column per species:

    Desc <tab> Family ID <tab> species_1 <tab> species_2 ...

Two filters are applied, following the CAFE5 documentation:

  1. Families that are absent from one of the two clades descending from the
     root are removed. CAFE5 assumes every family was present in the most
     recent common ancestor of all taxa, and families that fail this cannot be
     modelled sensibly. With Tribolium castaneum as the outgroup, this means a
     family must have at least one gene in T. castaneum and at least one gene
     in the ingroup.

  2. Families with a very large count in any single species are written to a
     separate file. Extreme families destabilise the birth-death rate estimate,
     so they are normally analysed in a second run with their own lambda.

The optional rename map is a two-column TSV, "old<tab>new", used to make the
count-table headers match the tip labels of the species tree exactly. CAFE5
matches species by name, and any mismatch is a fatal error.

Run from the repository root:
    conda activate oimp_05_gene_family_evolution
    python 05_gene_family_evolution/01_prepare_cafe_input.py
"""

import os
from types import SimpleNamespace


import csv
import sys

# --------------------------------- settings ----------------------------------
SETTINGS = SimpleNamespace(
    counts=os.environ.get(
        "COUNTS", "results/04_orthology/orthofinder_out/Orthogroups.GeneCount.tsv"),
    outdir=os.environ.get("OUTDIR", "results/05_gene_family_evolution"),
    outgroup=os.environ.get("OUTGROUP", "Tribolium_castaneum"),
    max_count=int(os.environ.get("MAX_COUNT", 30)),   # families with this many
                                                      # genes or more in any one
                                                      # species go to a separate
                                                      # file; CAFE cannot fit them
    rename=os.environ.get("RENAME") or None,          # optional two-column TSV
)
# -----------------------------------------------------------------------------


def read_counts(path):
    """Read OrthoFinder Orthogroups.GeneCount.tsv (or any equivalent table)."""
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter='\t'))
    header = rows[0]
    # OrthoFinder writes an id column first and a 'Total' column last
    id_col = 0
    species = [h for h in header[1:] if h.strip().lower() != 'total']
    idx = {h: i for i, h in enumerate(header)}
    out = []
    for r in rows[1:]:
        if not r or not r[id_col].strip():
            continue
        counts = {}
        for s in species:
            v = r[idx[s]].strip()
            counts[s] = int(v) if v.isdigit() else 0
        out.append((r[id_col].strip(), counts))
    return species, out


def main():
    args = SETTINGS

    species, families = read_counts(args.counts)

    if args.rename:
        mapping = {}
        with open(args.rename) as fh:
            for line in fh:
                if line.strip() and not line.startswith('#'):
                    old, new = line.rstrip('\n').split('\t')[:2]
                    mapping[old] = new
        species = [mapping.get(s, s) for s in species]
        families = [(og, {mapping.get(k, k): v for k, v in c.items()})
                    for og, c in families]

    if args.outgroup not in species:
        sys.exit(f'ERROR: outgroup {args.outgroup!r} is not a column.\n'
                 f'Columns are: {", ".join(species)}')
    ingroup = [s for s in species if s != args.outgroup]

    keep, large, dropped = [], [], 0
    for og, c in families:
        if c[args.outgroup] == 0 or all(c[s] == 0 for s in ingroup):
            dropped += 1
            continue
        (large if max(c.values()) >= args.max_count else keep).append((og, c))

    os.makedirs(args.outdir, exist_ok=True)

    def write(path, rows):
        with open(path, 'w') as fh:
            fh.write('Desc\tFamily ID\t' + '\t'.join(species) + '\n')
            for og, c in rows:
                fh.write('(null)\t' + og + '\t' +
                         '\t'.join(str(c[s]) for s in species) + '\n')
        return path

    f_main = write(os.path.join(args.outdir, 'cafe_gene_counts_filtered.tsv'), keep)
    f_large = write(os.path.join(args.outdir, 'cafe_gene_counts_large.tsv'), large)

    print(f'species                              : {len(species)}')
    print(f'orthogroups read                     : {len(families)}')
    print(f'  absent from one side of the root   : {dropped}  (removed)')
    print(f'  >= {args.max_count} genes in one species        : {len(large)}  '
          f'(separate run)')
    print(f'  retained for the main analysis     : {len(keep)}')
    print(f'  genes in retained families         : '
          f'{sum(sum(c.values()) for _, c in keep):,}')
    print()
    print('written:', f_main)
    print('written:', f_large)
    print()
    print('Species columns, in order (these must match the tree tip labels exactly):')
    for s in species:
        print('  ', s)


if __name__ == '__main__':
    main()
