#!/usr/bin/env bash
#
# Partitioned maximum-likelihood tree from the supermatrix (Fig. 2).
#
# -m MFP+MERGE runs ModelFinder and then merges partitions with similar models,
# which reduced the 593 gene partitions to 41 schemes. Support comes from 1,000
# ultrafast bootstrap replicates and 1,000 SH-aLRT tests.
#
# Run from the repository root:
#   conda activate oimp_04_orthology
#   bash 04_orthology/06_iqtree.sh
#
set -euo pipefail

# ------------------------------- settings ------------------------------------
PROJECT="${PROJECT:-results/04_orthology}"           # section output root
CONCAT="${CONCAT:-${PROJECT}/concat}"                # supermatrix, from step 05
IQTREE_DIR="${IQTREE_DIR:-${PROJECT}/iqtree}"        # output directory
# -----------------------------------------------------------------------------

mkdir -p "$IQTREE_DIR"
SUPERMATRIX="${CONCAT}/supermatrix_aa.fasta"
PARTITIONS="${CONCAT}/partitions_aa.txt"

command -v iqtree >/dev/null || { echo "Error: iqtree not found in PATH"; exit 1; }

# 1) Fit the partitioned model and search. IQ-TREE writes into the working
#    directory, so run it from there.
echo "Running IQ-TREE..."
cd "$IQTREE_DIR"

iqtree \
  -s "$SUPERMATRIX" \
  -p "$PARTITIONS" \
  -m MFP+MERGE \
  -bb 1000 -alrt 1000 \
  -nt AUTO

# 2) Report where the tree and its support values ended up.
echo "Done. Key outputs in $IQTREE_DIR:"
echo "  *.treefile   = Best ML tree"
echo "  *.ufboot     = Bootstrap trees"
echo "  *.iqtree     = Run report"
echo "  *.log        = Log file"
