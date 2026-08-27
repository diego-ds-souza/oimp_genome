#!/usr/bin/env bash
#
# Parse the curated GFF3 into the three tables the expression analyses need:
#
#   tx2gene.tsv                    transcript to gene, for tximport
#   gene_family_assignments.tsv    family membership per gene, from the Pfam,
#                                  MEROPS and CAZy annotations (Table S3)
#   families/<FAMILY>.txt          one gene list per family
#
# Family membership is read from the mRNA attributes of the annotation, so it
# reflects exactly what was deposited rather than a separate curation step.
#
# Run from the repository root:
#   conda activate oimp_07_transcriptomics
#   bash 07_transcriptomics/01_parse_gff3_annotations.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
GFF3="${GFF3:-results/03_annotation/curated/Oncideres_impluviata.clean_AGAT.renamed.products.fixed.renumbered.gff3}"
OUTDIR="${OUTDIR:-results/07_transcriptomics/ref}"        # tx2gene and gene names
FAMDIR="${FAMDIR:-results/07_transcriptomics/families}"   # one file per family
# -----------------------------------------------------------------------------

LOGDIR="results/07_transcriptomics/logs"
mkdir -p "$OUTDIR" "$FAMDIR" "$LOGDIR"

STAMP=$(date +"%Y%m%d_%H%M%S")
LOG="${LOGDIR}/01_parse_gff3_${STAMP}.log"

echo "=============================================================================" | tee "$LOG"
echo "Parse GFF3 Annotations" | tee -a "$LOG"
echo "=============================================================================" | tee -a "$LOG"
echo "Input GFF3: $GFF3" | tee -a "$LOG"
echo "Output dir: $OUTDIR" | tee -a "$LOG"
echo "Families:   $FAMDIR" | tee -a "$LOG"
echo "Started:    $(date)" | tee -a "$LOG"
echo "" | tee -a "$LOG"

# Verify input exists
if [[ ! -f "$GFF3" ]]; then
    echo "[ERROR] GFF3 file not found: $GFF3" | tee -a "$LOG"
    exit 1
fi

# 1) Create temporary python script (avoids heredoc issues with spaces in paths).
PYSCRIPT=$(mktemp)
trap "rm -f $PYSCRIPT" EXIT

cat > "$PYSCRIPT" << 'PYEOF'
import sys
import re
import os
from collections import defaultdict

gff3_path = sys.argv[1]
out_dir = sys.argv[2]
fam_dir = sys.argv[3]
log_path = sys.argv[4]

def log(msg):
    print(msg)
    with open(log_path, 'a') as f:
        f.write(msg + '\n')

FAMILY_PFAM = {
    'OBP':   ['PF01395'],
    'CSP':   ['PF03392'],
    'OR':    ['PF02949', 'PF13853'],
    'GR':    ['PF08395', 'PF06151'],
    'IR':    ['PF00060', 'PF10613'],
    'SNMP':  ['PF01130'],
    'P450':  ['PF00067'],
    'GST':   ['PF00043', 'PF02798', 'PF13417', 'PF13410', 'PF14497'],
    'UGT':   ['PF00201'],
    'CCE':   ['PF00135'],
    'ABC':   ['PF00005', 'PF00664', 'PF06472', 'PF00004'],
    'SerPro': ['PF00089', 'PF00082', 'PF05577'],
    'CysPro': ['PF00112', 'PF00656'],
    'AspPro': ['PF00026', 'PF07966'],
    'MetPro': ['PF01546', 'PF00557'],
    'FAS':   ['PF00109', 'PF02801', 'PF00698'],
    'ELO':   ['PF01151'],
    'Desat': ['PF00487'],
    'ACSL':  ['PF00501', 'PF13193'],
    'FPPS':  ['PF00348'],
    'IspS':  ['PF01397'],
}

MEROPS_SERINE = ['S01', 'S08', 'S09', 'S10', 'S28', 'S33']

CAZY_FAMILIES = {
    'GH':  r'GH\d+',
    'GT':  r'GT\d+',
    'PL':  r'PL\d+',
    'CE':  r'CE\d+',
    'CBM': r'CBM\d+',
    'AA':  r'AA\d+',
}

log("[Step 1] Parsing GFF3 file...")

tx2gene = {}
gene_annotations = defaultdict(lambda: {
    'pfam': set(), 'interpro': set(), 'go': set(),
    'merops': set(), 'cazy': set(), 'eggnog': set(),
    'product': '', 'name': ''
})

n_genes = 0
n_mrna = 0

with open(gff3_path, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        cols = line.strip().split('\t')
        if len(cols) < 9:
            continue
        feature_type = cols[2]
        attributes = cols[8]
        
        if feature_type == 'gene':
            n_genes += 1
            name_match = re.search(r'Name=([^;]+)', attributes)
            id_match = re.search(r'ID=([^;]+)', attributes)
            if id_match:
                gid = id_match.group(1)
                if name_match:
                    gene_annotations[gid]['name'] = name_match.group(1)
            continue
        
        if feature_type == 'mRNA':
            n_mrna += 1
            tid_match = re.search(r'ID=([^;]+)', attributes)
            parent_match = re.search(r'Parent=([^;]+)', attributes)
            if not tid_match or not parent_match:
                continue
            tid = tid_match.group(1)
            gid = parent_match.group(1)
            tx2gene[tid] = gid
            
            dbxref_match = re.search(r'Dbxref=([^;]+)', attributes)
            if dbxref_match:
                for ref in dbxref_match.group(1).split(','):
                    if ref.startswith('PFAM:'):
                        gene_annotations[gid]['pfam'].add(ref.replace('PFAM:', ''))
                    elif ref.startswith('InterPro:'):
                        gene_annotations[gid]['interpro'].add(ref.replace('InterPro:', ''))
            
            go_match = re.search(r'Ontology_term=([^;]+)', attributes)
            if go_match:
                for go in go_match.group(1).split(','):
                    if go.startswith('GO:'):
                        gene_annotations[gid]['go'].add(go)
            
            note_match = re.search(r'note=([^;]+)', attributes)
            if note_match:
                for item in note_match.group(1).split(','):
                    if item.startswith('EggNog:'):
                        gene_annotations[gid]['eggnog'].add(item.replace('EggNog:', ''))
                    elif item.startswith('MEROPS:'):
                        gene_annotations[gid]['merops'].add(item.replace('MEROPS:', ''))
                    elif item.startswith('CAZy:'):
                        gene_annotations[gid]['cazy'].add(item.replace('CAZy:', ''))
            
            product_match = re.search(r'product=([^;]+)', attributes)
            if product_match:
                gene_annotations[gid]['product'] = product_match.group(1)

log(f"  Genes found: {n_genes}")
log(f"  mRNAs found: {n_mrna}")
log(f"  Unique transcripts: {len(tx2gene)}")

log("\n[Step 2] Writing tx2gene.tsv...")
tx2gene_path = os.path.join(out_dir, 'tx2gene.tsv')
with open(tx2gene_path, 'w') as f:
    f.write("transcript_id\tgene_id\n")
    for tid, gid in sorted(tx2gene.items()):
        f.write(f"{tid}\t{gid}\n")
log(f"  Wrote: {tx2gene_path} ({len(tx2gene)} transcripts)")

log("\n[Step 3] Assigning genes to families...")
gene_families = defaultdict(set)

for gid, annot in gene_annotations.items():
    pfams = annot['pfam']
    merops = annot['merops']
    cazy = annot['cazy']
    
    for family, pfam_list in FAMILY_PFAM.items():
        if pfams & set(pfam_list):
            gene_families[gid].add(family)
    
    for m in merops:
        for sp_fam in MEROPS_SERINE:
            if sp_fam in m:
                gene_families[gid].add('SerPro')
                break
    
    for cazy_entry in cazy:
        for cazy_name, pattern in CAZY_FAMILIES.items():
            if re.search(pattern, cazy_entry):
                gene_families[gid].add('CAZyme')
                gene_families[gid].add(f'CAZy_{cazy_name}')

log("\n[Step 3b] Resolving dual-assigned genes...")

for gid in gene_families:
    fams = gene_families[gid]
    if 'UGT' in fams:
        fams -= {'CAZyme', 'CAZy_GT'}
    if 'SerPro' in fams:
        fams -= {'CAZyme', 'CAZy_CBM', 'CAZy_GH'}

# Count and report resolved conflicts
n_resolved = 0
for gid in gene_families:
    fams = gene_families[gid]
    if 'UGT' in fams and fams & {'CAZyme', 'CAZy_GT'}:
        n_resolved += 1
    if 'SerPro' in fams and fams & {'CAZyme', 'CAZy_CBM', 'CAZy_GH'}:
        n_resolved += 1

log(f"  Dual assignments resolved")

log("\n[Step 4] Writing gene family assignments...")
assign_path = os.path.join(fam_dir, 'gene_family_assignments.tsv')
with open(assign_path, 'w') as f:
    f.write("gene_id\tfamilies\tpfam\tmerops\tcazy\tproduct\n")
    for gid in sorted(gene_annotations.keys()):
        annot = gene_annotations[gid]
        families = ','.join(sorted(gene_families.get(gid, set()))) or 'NA'
        pfams = ','.join(sorted(annot['pfam'])) or 'NA'
        merops = ','.join(sorted(annot['merops'])) or 'NA'
        cazy = ','.join(sorted(annot['cazy'])) or 'NA'
        product = annot['product'] or 'NA'
        f.write(f"{gid}\t{families}\t{pfams}\t{merops}\t{cazy}\t{product}\n")
log(f"  Wrote: {assign_path}")

all_families = set()
for fams in gene_families.values():
    all_families.update(fams)

family_counts = {}
for fam in sorted(all_families):
    genes = [gid for gid, fams in gene_families.items() if fam in fams]
    family_counts[fam] = len(genes)
    fam_path = os.path.join(fam_dir, f'genes_{fam}.txt')
    with open(fam_path, 'w') as f:
        for g in sorted(genes):
            f.write(f"{g}\n")

log("\n[Step 5] Family counts:")
log("-" * 40)

chemosensory = ['OBP', 'CSP', 'OR', 'GR', 'IR', 'SNMP']
detox = ['P450', 'GST', 'UGT', 'CCE', 'ABC']
digestion = ['SerPro', 'CysPro', 'AspPro', 'MetPro', 'CAZyme']
lipid_metabolism = ['FAS', 'ELO', 'Desat', 'ACSL']
terpenoid = ['FPPS', 'IspS']
cazy_specific = [f for f in all_families if f.startswith('CAZy_')]

log("\nChemosensory:")
for fam in chemosensory:
    if fam in family_counts:
        log(f"  {fam:8s}: {family_counts[fam]:4d} genes")

log("\nDetoxification:")
for fam in detox:
    if fam in family_counts:
        log(f"  {fam:8s}: {family_counts[fam]:4d} genes")

log("\nDigestion (proteases + CAZymes):")
for fam in digestion:
    if fam in family_counts and not fam.startswith('CAZy_'):
        log(f"  {fam:8s}: {family_counts[fam]:4d} genes")

log("\nLipid metabolism:")
for fam in lipid_metabolism:
    if fam in family_counts:
        log(f"  {fam:8s}: {family_counts[fam]:4d} genes")

log("\nTerpenoid/Isoprenoid metabolism:")
for fam in terpenoid:
    if fam in family_counts:
        log(f"  {fam:8s}: {family_counts[fam]:4d} genes")

log("\nCAZyme classes:")
for fam in sorted(cazy_specific):
    if fam in family_counts:
        log(f"  {fam:10s}: {family_counts[fam]:4d} genes")

total_assigned = len([g for g in gene_families if gene_families[g]])
log(f"\n{'-' * 40}")
log(f"Total genes in GFF3:         {n_genes}")
log(f"Genes assigned to >=1 family: {total_assigned}")

log("\n[Step 6] Writing gene names map...")
names_path = os.path.join(out_dir, 'gene_names.tsv')
with open(names_path, 'w') as f:
    f.write("gene_id\tname\tproduct\n")
    for gid in sorted(gene_annotations.keys()):
        annot = gene_annotations[gid]
        name = annot['name'] or gid
        product = annot['product'] or 'hypothetical protein'
        f.write(f"{gid}\t{name}\t{product}\n")
log(f"  Wrote: {names_path}")

log("\n[Step 7] Writing gene2go mapping...")
g2go_path = os.path.join(out_dir, 'gene2go.tsv')
n_go_pairs = 0
with open(g2go_path, 'w') as f:
    for gid in sorted(gene_annotations.keys()):
        for go in sorted(gene_annotations[gid]['go']):
            f.write(f"{gid}\t{go}\n")
            n_go_pairs += 1

genes_with_go = sum(1 for g in gene_annotations.values() if g['go'])
log(f"  Wrote: {g2go_path}")
log(f"  Genes with GO annotations: {genes_with_go}")
log(f"  Total gene-GO pairs: {n_go_pairs}")

log("\n" + "=" * 70)
log("DONE!")
log("=" * 70)
PYEOF

# 2) Run the python script with properly quoted arguments.
python3 "$PYSCRIPT" "$GFF3" "$OUTDIR" "$FAMDIR" "$LOG"

echo "" | tee -a "$LOG"
echo "Output files:" | tee -a "$LOG"
echo "  - ${OUTDIR}/tx2gene.tsv" | tee -a "$LOG"
echo "  - ${OUTDIR}/gene_names.tsv" | tee -a "$LOG"
echo "  - ${OUTDIR}/gene2go.tsv" | tee -a "$LOG"
echo "  - ${FAMDIR}/gene_family_assignments.tsv" | tee -a "$LOG"
echo "  - ${FAMDIR}/genes_*.txt (per-family gene lists)" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "Completed: $(date)" | tee -a "$LOG"
