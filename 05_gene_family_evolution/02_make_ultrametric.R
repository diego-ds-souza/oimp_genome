#!/usr/bin/env Rscript
#
# Convert the IQ-TREE species tree into the chronogram CAFE5 requires.
#
# The IQ-TREE result is a phylogram: branch lengths are substitutions per site,
# not time, and the tips are not equidistant from the root. CAFE5 needs a
# rooted, binary, ultrametric tree.
#
# Three things happen here:
#   the tree is rooted on the outgroup and node labels are dropped, since
#   bootstrap values confuse some downstream parsers; ape::chronos() applies
#   penalized-likelihood rate smoothing to produce a chronogram; and the result
#   is scaled so the root sits at the calibration age.
#
# chronos() can leave a root child at the same age as the root, giving a branch
# of ~1e-14 that CAFE5 cannot evaluate. Step 03 repairs that; step 04 checks it.
#
# Run from the repository root:
#   Rscript 05_gene_family_evolution/02_make_ultrametric.R

# ------------------------------- settings ------------------------------------
RES      <- Sys.getenv("RES", "results/05_gene_family_evolution")
infile   <- file.path("results/04_orthology/concat",
                      "partitions_aa.txt_rooted_v2.treefile")
outfile  <- file.path(RES, "species_tree_ultrametric.nwk")
outgroup <- "Tribolium_castaneum"
root_age <- 227          # Ma; Phytophaga crown age used to calibrate the root
# -----------------------------------------------------------------------------

suppressPackageStartupMessages(library(ape))
dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)

tr <- read.tree(infile)
message("Tips read: ", Ntip(tr))

if (!outgroup %in% tr$tip.label) {
  stop("Outgroup '", outgroup, "' is not a tip label.\nTips are:\n  ",
       paste(tr$tip.label, collapse = "\n  "))
}

# ---- root, resolve polytomies, strip node labels -----------------------------
tr <- root(tr, outgroup = outgroup, resolve.root = TRUE)
if (!is.binary(tr)) {
  message("Tree was not binary; resolving polytomies with multi2di()")
  tr <- multi2di(tr)
}
tr$node.label <- NULL

# ---- rate smoothing ----------------------------------------------------------
message("Running chronos() ...")
chr <- chronos(tr, lambda = 1, model = "correlated",
               control = chronos.control(nb.rate.cat = 1))
class(chr) <- "phylo"

# ---- scale the root ----------------------------------------------------------
depth <- max(node.depth.edgelength(chr))
chr$edge.length <- chr$edge.length * (root_age / depth)

stopifnot(is.rooted(chr), is.binary(chr))
if (!is.ultrametric(chr, tol = 1e-6)) {
  stop("chronos() did not return an ultrametric tree; inspect the input")
}

write.tree(chr, file = outfile)

message("")
message("Root depth set to : ", root_age,
        if (length(args) >= 4) " Ma" else " (relative units)")
message("Ultrametric       : ", is.ultrametric(chr))
message("Binary and rooted : ", is.binary(chr), " / ", is.rooted(chr))
message("Written           : ", outfile)
message("")
message("Tip labels, in order -- these must match the count-table headers exactly:")
for (t in sort(chr$tip.label)) message("  ", t)
