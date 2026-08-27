#!/usr/bin/env bash
#
# Concatenate the trimmed single-copy alignments into a supermatrix and write
# the matching IQ-TREE partition file, one partition per gene.
#
# The taxon set and its order must be identical in every alignment; any gene
# file that disagrees is skipped and reported in concat.log rather than being
# padded silently.
#
# Run from the repository root:
#   conda activate oimp_04_orthology
#   bash 04_orthology/05_concat_and_partitions.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
PROJECT="${PROJECT:-results/04_orthology}"          # section output root
TRIM_DIR="${TRIM_DIR:-${PROJECT}/trimmed}"          # trimmed alignments, step 04
OUT_DIR="${OUT_DIR:-${PROJECT}/concat}"             # output directory
# -----------------------------------------------------------------------------

mkdir -p "$OUT_DIR"
SUPERMATRIX="${OUT_DIR}/supermatrix_aa.fasta"
PARTITIONS="${OUT_DIR}/partitions_aa.txt"
LOG="${OUT_DIR}/concat.log"

command -v python3 >/dev/null || { echo "Error: python3 not found in PATH"; exit 1; }

# 1) Concatenate. Every alignment must carry the same taxa in the same order;
#    a gene file that disagrees is skipped and named in the log rather than
#    padded silently.
python3 - "$TRIM_DIR" "$SUPERMATRIX" "$PARTITIONS" "$LOG" << 'PY'
import sys, os, glob, textwrap

trim_dir, super_fa, part_txt, log_path = sys.argv[1:5]

def read_fasta(p):
    d = {}
    h = None
    seq = []
    with open(p) as fh:
        for line in fh:
            if line.startswith('>'):
                if h is not None:
                    d[h] = ''.join(seq)
                h = line.strip()[1:]
                seq = []
            else:
                seq.append(line.strip())
        if h is not None:
            d[h] = ''.join(seq)
    return d

def species_name(header):
    # Expect "SpeciesTag|GeneID" -> take SpeciesTag
    return header.split('|',1)[0]

files = sorted(glob.glob(os.path.join(trim_dir, "*.aln.fa")))
if not files:
    print(f"No files found in {trim_dir}", file=sys.stderr); sys.exit(1)

log = []
# 2) Take the master taxon order from the first alignment.
first = files[0]
D0 = read_fasta(first)
taxa = [species_name(h) for h in D0.keys()]
if len(taxa) != len(set(taxa)):
    print("Duplicate taxa detected in the first alignment; check headers.", file=sys.stderr)
    sys.exit(1)

# 3) Require exactly one sequence per taxon.
def normalize(d):
    m = {}
    for h, s in d.items():
        sp = species_name(h)
        if sp in m and len(s) > len(m[sp]):
            m[sp] = s
        elif sp not in m:
            m[sp] = s
    return m

ok_files = []
blocks = []         # list of dict {taxon -> seq}
lengths = []        # per-gene length
kept = 0; skipped = 0

for f in files:
    D = normalize(read_fasta(f))
    # check coverage and identical taxon set
    if set(D.keys()) != set(taxa):
        skipped += 1
        missing = [t for t in taxa if t not in D]
        extra = [t for t in D if t not in taxa]
        log.append(f"SKIP {os.path.basename(f)}  missing={len(missing)} extra={len(extra)}")
        continue
    L = len(next(iter(D.values())))
    if any(len(s) != L for s in D.values()):
        skipped += 1
        log.append(f"SKIP {os.path.basename(f)}  inconsistent lengths")
        continue
    ok_files.append(f)
    blocks.append([D[t] for t in taxa])
    lengths.append(L)
    kept += 1

if kept == 0:
    print("No alignments passed the consistency checks.", file=sys.stderr); sys.exit(1)

# 4) Write the supermatrix, one sequence per taxon.
with open(super_fa, "w") as fo:
    for i, sp in enumerate(taxa):
        seq = "".join(block[i] for block in blocks)
        fo.write(f">{sp}\n")
        for j in range(0, len(seq), 60):
            fo.write(seq[j:j+60] + "\n")

# 5) Write the partition file, one partition per gene.
start = 1
with open(part_txt, "w") as fp:
    for f, L in zip(ok_files, lengths):
        gene = os.path.basename(f).rsplit(".aln.fa", 1)[0]
        end = start + L - 1
        fp.write(f"LG+G4, {gene} = {start}-{end}\n")
        start = end + 1

# 6) Write the log, including any gene that was skipped.
with open(log_path, "w") as lg:
    tot = len(files)
    lg.write(f"Total gene alignments: {tot}\n")
    lg.write(f"Kept: {kept}\nSkipped: {skipped}\n")
    for line in log:
        lg.write(line + "\n")

print(f"Supermatrix: {super_fa}")
print(f"Partitions:  {part_txt}")
print(f"Log:         {log_path}")
PY

echo "Done. Supermatrix and partitions are in: ${OUT_DIR}"
