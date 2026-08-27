#!/usr/bin/env Rscript
#
# Diagnose a tree that CAFE5 refuses to use.
#
# CAFE5 reports -lnL = inf on every lambda when a branch length is numerically
# zero. Such a branch is positive, so a naive check misses it; this prints the
# smallest branch, the tip-depth spread and the counts that reveal it.
#
# Run it on the repaired tree from step 03, and on the unrepaired tree from
# step 02 for comparison.
#
# Run from the repository root:
#   Rscript 05_gene_family_evolution/04_check_tree.R

# ------------------------------- settings ------------------------------------
RES  <- Sys.getenv("RES", "results/05_gene_family_evolution")
TREE <- Sys.getenv("TREE", file.path(RES, "species_tree_ultrametric_fixed.nwk"))
# -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(ape))

tr <- read.tree(TREE)
bl <- tr$edge.length
cat("tips              :", Ntip(tr), "\n")
cat("rooted            :", is.rooted(tr), "\n")
cat("binary            :", is.binary(tr), "\n")
cat("ultrametric 1e-6  :", is.ultrametric(tr, tol=1e-6), "\n")
cat("ultrametric 1e-8  :", is.ultrametric(tr, tol=1e-8), "\n")
cat("node labels       :", if (is.null(tr$node.label)) "none" else "PRESENT (remove)", "\n")
cat("branches          :", length(bl), "\n")
cat("  min             :", min(bl), "\n")
cat("  max             :", max(bl), "\n")
cat("  zero-length     :", sum(bl == 0), "\n")
cat("  negative        :", sum(bl < 0), "\n")
cat("  < 1e-6          :", sum(bl < 1e-6), "\n")
d <- node.depth.edgelength(tr)[1:Ntip(tr)]
cat("root-to-tip depth : min", min(d), " max", max(d),
    " spread", max(d)-min(d), "\n")
if (sum(bl <= 0) > 0 || (max(d)-min(d)) > 1e-6)
  cat("\n*** This is very likely why CAFE5 returns -lnL = inf ***\n")
