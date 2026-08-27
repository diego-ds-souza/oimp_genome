#!/usr/bin/env python3
"""
Draw the CAFE expansion/contraction column proposed for Figure 2, in two
encodings, so the choice can be made from the real numbers:

  <prefix>_bars.pdf / .png    diverging bar on a shared scale   (recommended)
  <prefix>_pies.pdf / .png    pie per species                   (the convention
                                                                 in many papers)
  <prefix>_values.tsv         the numbers, in tree order, for grafting the
                              column onto the ggtree figure in R

Both figures show the species name, the mark, and the "+expanded / -contracted"
values, one row per tip, in the same order as the tree.

Why two encodings: a pie shows only the ratio, so a branch with 40 expansions
and one with 1,340 look equally important, and ratios in the 60-90% range are
hard to rank as angles. The diverging bar puts every row on one shared scale and
keeps magnitude visible. The pie version is provided because it is the familiar
convention; compare them at print size before deciding.

Usage
  python3 05_plot_cafe_column.py \
      --clade results/base_error/Base_clade_results.txt \
      --tree  results/partitions_aa_ultrametric_227_fixed.nwk \
      --outdir figures --prefix cafe_column

Options worth knowing
  --order-file     text file, one species per line, top row first. This is the
                   reliable way to match Figure 2: R's plotting order depends on
                   the tree's internal edge order after ladderize(), which the
                   Newick string only approximates. For the R version of the
                   figure you do not need it - join the values TSV to the tree
                   data by species name and let ggtree supply the order.
  --no-reverse     take the Newick tip order as-is instead of reversing it.
  --row-height     inches per species row; raise it to match Figure 2's spacing
  --focal          species drawn in bold (default Oncideres_impluviata)

Run from the repository root:
    conda activate oimp_05_gene_family_evolution
    python 05_gene_family_evolution/06_plot_cafe_column.py
"""

import os
from types import SimpleNamespace


import csv
import os
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --------------------------------- settings ----------------------------------
RES = os.environ.get("RES", "results/05_gene_family_evolution")
SETTINGS = SimpleNamespace(
    clade=os.environ.get("CLADE", os.path.join(RES, "base_error",
                                               "Base_clade_results.txt")),
    tree=os.environ.get("TREE", os.path.join(RES,
                                             "species_tree_ultrametric_fixed.nwk")),
    outdir=os.environ.get("OUTDIR", os.path.join(RES, "plots")),
    prefix=os.environ.get("PREFIX", "cafe_column"),
    focal=os.environ.get("FOCAL", "Oncideres_impluviata"),
    row_height=float(os.environ.get("ROW_HEIGHT", 0.33)),
    width=float(os.environ.get("WIDTH", 7.0)),
    no_reverse=bool(os.environ.get("NO_REVERSE")),
    order_file=os.environ.get("ORDER_FILE") or None,   # tip order, if not the tree's
)
# -----------------------------------------------------------------------------

# Palette for expansions and contractions.
EXP = '#3C8DBC'      # expansions
CON = '#E8A33D'      # contractions
INK = '#111111'
MUTED = '#5A5A5A'
RULE = '#999999'

# Display names. The pipeline uses whatever the FASTA files were called; the
# manuscript uses these. Edit freely - anything not listed is just printed with
# underscores turned into spaces.
DISPLAY = {
    'Anthonomus_grandis': 'Anthonomus grandis grandis',
    'Diabrotica_virgifera': 'Diabrotica virgifera virgifera',
    'Psylliodes_chrysocephalus': 'Psylliodes chrysocephala',
}


def tip_order(tree_path, reverse=True):
    """Tip labels in the order they appear in the Newick file."""
    nwk = open(tree_path).read()
    tips = re.findall(r'[(,]\s*([A-Za-z0-9_.\-]+)\s*:', nwk)
    seen, out = set(), []
    for t in tips:
        if t not in seen:
            seen.add(t)
            out.append(t)
    return out[::-1] if reverse else out


def read_clade(path):
    """CAFE *_clade_results.txt -> {name: (increase, decrease)}, tips only."""
    tips, internal = {}, {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 3 or not f[0].strip():
                continue
            m = re.match(r'^(.*?)<(\d+)>$', f[0].strip())
            if not m:
                continue
            name, nid = m.group(1), int(m.group(2))
            rec = (int(f[1]), int(f[2]), nid)
            (internal if name == '' else tips)[name or f'node_{nid}'] = rec
    return tips, internal


def label(name):
    return DISPLAY.get(name, name.replace('_', ' '))


# Layout shared by both encodings.
def new_figure(rows, row_height, width):
    n = len(rows)
    fig, ax = plt.subplots(figsize=(width, max(2.0, n * row_height + 1.0)))
    ax.set_xlim(0, 10)
    ax.set_ylim(-1.3, n)
    ax.axis('off')
    return fig, ax, n


def draw_names_and_values(ax, rows, n, focal, x_name=0.0, x_val=9.95):
    for i, (name, e, c, _) in enumerate(rows):
        y = n - 1 - i
        bold = name == focal
        ax.text(x_name, y, label(name), fontsize=8.6, style='italic',
                color=INK, va='center', ha='left',
                fontweight='bold' if bold else 'normal')
        ax.text(x_val, y, f'+{e:,} / −{c:,}', fontsize=8.2, color=INK,
                va='center', ha='right',
                fontweight='bold' if bold else 'normal')


def plot_bars(rows, focal, row_height, width, out):
    fig, ax, n = new_figure(rows, row_height, width)
    draw_names_and_values(ax, rows, n, focal)
    mx = max(max(e, c) for _, e, c, _ in rows)
    zero, half = 6.4, 1.85
    for i, (name, e, c, _) in enumerate(rows):
        y = n - 1 - i
        ax.barh(y, half * e / mx, left=zero, height=0.56, color=EXP, zorder=3)
        ax.barh(y, -half * c / mx, left=zero, height=0.56, color=CON, zorder=3)
    ax.plot([zero, zero], [-0.5, n - 0.55], color=RULE, lw=0.8, zorder=4)
    for s, lab in ((-half, f'{mx:,}'), (0, '0'), (half, f'{mx:,}')):
        ax.text(zero + s, -0.85, lab, fontsize=7.2, color=MUTED,
                ha='center', va='center')
    ax.text(zero - half / 2, -1.22, 'contracted', fontsize=7.6, color=CON,
            ha='center', va='center')
    ax.text(zero + half / 2, -1.22, 'expanded', fontsize=7.6, color=EXP,
            ha='center', va='center')
    save(fig, out)


def plot_pies(rows, focal, row_height, width, out):
    fig, ax, n = new_figure(rows, row_height, width)
    draw_names_and_values(ax, rows, n, focal)
    cx, r = 6.4, 0.36
    for i, (name, e, c, _) in enumerate(rows):
        y = n - 1 - i
        frac = e / (e + c) if (e + c) else 0.0
        pa = ax.inset_axes([cx - r, y - r, 2 * r, 2 * r],
                           transform=ax.transData)
        pa.set_aspect('equal')
        pa.pie([frac, 1 - frac], colors=[EXP, CON], startangle=90,
               counterclock=False,
               wedgeprops=dict(edgecolor='white', linewidth=0.7))
        pa.set_axis_off()
    ax.text(cx, -1.0, 'blue: expanded   orange: contracted', fontsize=7.6,
            color=MUTED, ha='center', va='center')
    ax.text(cx, -1.28, 'area is a proportion; rows are not comparable',
            fontsize=7.0, color=MUTED, ha='center', va='center')
    save(fig, out)


def save(fig, stem):
    fig.tight_layout()
    for ext in ('pdf', 'png'):
        p = f'{stem}.{ext}'
        fig.savefig(p, dpi=400, facecolor='white')
        print('wrote', p)
    plt.close(fig)


def main():
    args = SETTINGS
    os.makedirs(args.outdir, exist_ok=True)

    if args.order_file:
        order = [l.strip() for l in open(args.order_file) if l.strip()]
        print(f'row order taken from {args.order_file}')
    else:
        order = tip_order(args.tree, reverse=not args.no_reverse)
    tips, internal = read_clade(args.clade)

    missing = [t for t in order if t not in tips]
    extra = [t for t in tips if t not in order]
    if missing or extra:
        raise SystemExit(
            'FATAL: tree tips and clade-result tips disagree.\n'
            + (f'  in tree, not in clade file: {", ".join(missing)}\n' if missing else '')
            + (f'  in clade file, not in tree: {", ".join(extra)}\n' if extra else ''))

    rows = [(t, tips[t][0], tips[t][1], tips[t][2]) for t in order]
    if args.order_file:
        print(f'{len(rows)} tips, order from --order-file')
    else:
        print(f'{len(rows)} tips, in Newick order'
              f'{"" if args.no_reverse else ", reversed"}'
              '  (check this against Figure 2; if it is upside-down pass '
              '--no-reverse, if it differs in detail use --order-file)')

    stem = os.path.join(args.outdir, args.prefix)
    plot_bars(rows, args.focal, args.row_height, args.width, stem + '_bars')
    plot_pies(rows, args.focal, args.row_height, args.width, stem + '_pies')

    # values, in the same order, for the R version of the figure
    tsv = stem + '_values.tsv'
    with open(tsv, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['row', 'species', 'display_name', 'cafe_node',
                    'expanded', 'contracted', 'total', 'percent_expanded'])
        for i, (name, e, c, nid) in enumerate(rows, 1):
            tot = e + c
            w.writerow([i, name, label(name), nid, e, c, tot,
                        f'{100 * e / tot:.1f}' if tot else 'NA'])
    print('wrote', tsv)

    if internal:
        tsv2 = stem + '_internal_nodes.tsv'
        with open(tsv2, 'w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t')
            w.writerow(['cafe_node', 'expanded', 'contracted'])
            for k, (e, c, nid) in sorted(internal.items(), key=lambda x: x[1][2]):
                w.writerow([nid, e, c])
        print('wrote', tsv2, f'({len(internal)} internal branches)')

    foc = next((r for r in rows if r[0] == args.focal), None)
    if foc:
        print(f'\n{label(foc[0])}: +{foc[1]:,} / −{foc[2]:,} '
              f'({100 * foc[1] / (foc[1] + foc[2]):.1f}% expanded)')


if __name__ == '__main__':
    main()