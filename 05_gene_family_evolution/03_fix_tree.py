#!/usr/bin/env python3
"""
Repair an ultrametric tree that CAFE5 rejects with -lnL = inf on every lambda.

The cause is a branch whose length is numerically zero. `ape::root(...,
resolve.root = TRUE)` splits the root edge and can leave one root child at the
same age as the root; `multi2di()` does the same when resolving a polytomy.
After scaling the tree to a root depth of 250 the residue shows up as a branch
of ~1e-14, which is below double-precision resolution at that scale. Birth-death
transition probabilities contain terms that reduce to 0/0 as t -> 0, so the
matrix picks up a NaN and every family likelihood becomes inf, whatever lambda
is tried. The branch is not zero and not negative, so simple checks miss it.

The repair works in node-age space, so ultrametricity is preserved exactly:

  1. tips are at age 0, the root at the tree depth
  2. walking root-ward to tip-ward, any node whose parent branch is shorter than
     --min-branch is made younger, lengthening that branch
  3. every tip is then snapped to exactly age 0, removing residual root-to-tip
     jitter from the rate-smoothing step

Only internal node ages move, by a fraction of the tree depth. Topology is
untouched and no biological conclusion can depend on it.

Usage:
    python 09_b_fix_tree.py in.nwk out.nwk [--min-branch 0.25]

    python3 02_b_fix_tree.py \
  results/partitions_aa_ultrametric_227.nwk \
  results/partitions_aa_ultrametric_227_fixed.nwk \
  --min-branch 0.25

  
With no --min-branch the floor defaults to 0.1% of the tree depth.

Run from the repository root:
    conda activate oimp_05_gene_family_evolution
    python 05_gene_family_evolution/03_fix_tree.py
"""

import os
from types import SimpleNamespace


import re
import sys

# --------------------------------- settings ----------------------------------
RES = os.environ.get("RES", "results/05_gene_family_evolution")
SETTINGS = SimpleNamespace(
    infile=os.environ.get("INFILE", os.path.join(RES, "species_tree_ultrametric.nwk")),
    outfile=os.environ.get("OUTFILE",
                           os.path.join(RES, "species_tree_ultrametric_fixed.nwk")),
    min_branch=float(os.environ["MIN_BRANCH"]) if os.environ.get("MIN_BRANCH")
    else 0.25,        # Ma; branches shorter than this are lifted to it
)
# -----------------------------------------------------------------------------


# ----------------------------------------------------------------- parsing
class Node:
    __slots__ = ('name', 'length', 'children', 'parent', 'age')

    def __init__(self):
        self.name = ''
        self.length = 0.0
        self.children = []
        self.parent = None
        self.age = 0.0


def parse_newick(text):
    text = text.strip()
    if not text.endswith(';'):
        raise SystemExit('not a Newick string (no terminating ";")')
    tokens = re.findall(r'[(),;]|[^(),;]+', text)
    root, cur, i = Node(), None, 0
    stack = []
    cur = root
    for tok in tokens:
        if tok == '(':
            child = Node()
            child.parent = cur
            cur.children.append(child)
            stack.append(cur)
            cur = child
        elif tok == ',':
            parent = stack[-1]
            sib = Node()
            sib.parent = parent
            parent.children.append(sib)
            cur = sib
        elif tok == ')':
            cur = stack.pop()
        elif tok == ';':
            break
        else:
            label = tok.strip()
            if ':' in label:
                nm, _, ln = label.partition(':')
                cur.name = nm.strip()
                cur.length = float(ln)
            else:
                cur.name = label
    return root


def walk(node):
    yield node
    for c in node.children:
        yield from walk(c)


def format_newick(node, prec=12):
    if node.children:
        inner = ','.join(format_newick(c, prec) for c in node.children)
        s = f'({inner}){node.name}'
    else:
        s = node.name
    if node.parent is not None:
        s += f':{node.length:.{prec}f}'
    return s


# ------------------------------------------------------------------ repair
def depths(root):
    """root-to-node distance for every node"""
    d = {root: 0.0}
    for n in walk(root):
        for c in n.children:
            d[c] = d[n] + c.length
    return d


def main():
    args = SETTINGS

    root = parse_newick(open(args.infile).read())
    nodes = list(walk(root))
    tips = [n for n in nodes if not n.children]
    d = depths(root)
    depth = max(d[t] for t in tips)
    floor = args.min_branch if args.min_branch is not None else depth * 1e-3

    lens = [n.length for n in nodes if n.parent is not None]
    print(f'read              : {args.infile}')
    print(f'  tips            : {len(tips)}')
    print(f'  branches        : {len(lens)}')
    print(f'  tree depth      : {depth:g}')
    print(f'  min branch      : {min(lens):g}')
    print(f'  max branch      : {max(lens):g}')
    print(f'  root-tip spread : {max(d[t] for t in tips) - min(d[t] for t in tips):g}')
    print(f'  branches < {floor:g} : {sum(1 for L in lens if L < floor)}')

    # ages: tips at 0, root at depth
    for n in nodes:
        n.age = depth - d[n]

    # 1. lengthen short branches by making the child younger, root-ward first
    order = sorted((n for n in nodes if n.parent is not None),
                   key=lambda n: -n.age)
    moved = 0
    for n in order:
        if n.parent.age - n.age < floor:
            new_age = n.parent.age - floor
            if not n.children:
                continue                     # tips are pinned at 0, handled below
            n.age = new_age
            moved += 1
    # 2. pin tips to exactly 0
    for t in tips:
        t.age = 0.0

    # 3. recompute branch lengths from ages
    bad = 0
    for n in nodes:
        if n.parent is not None:
            n.length = n.parent.age - n.age
            if n.length < 0:
                bad += 1
    if bad:
        sys.exit(f'ERROR: {bad} negative branches after repair; '
                 f'lower --min-branch')

    d2 = depths(root)
    lens2 = [n.length for n in nodes if n.parent is not None]
    spread = max(d2[t] for t in tips) - min(d2[t] for t in tips)
    print()
    print(f'repaired          : {args.outfile}')
    print(f'  nodes moved     : {moved}')
    print(f'  min branch      : {min(lens2):g}')
    print(f'  max branch      : {max(lens2):g}')
    print(f'  tree depth      : {max(d2[t] for t in tips):g}')
    print(f'  root-tip spread : {spread:g}')
    if min(lens2) < floor * 0.999 or spread > 1e-9:
        sys.exit('ERROR: repair did not converge; inspect the tree')
    print('  OK: no branch below the floor, tips exactly level')

    with open(args.outfile, 'w') as fh:
        fh.write(format_newick(root) + ';\n')


if __name__ == '__main__':
    main()