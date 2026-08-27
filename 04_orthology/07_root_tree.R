#!/usr/bin/env Rscript
#
# Root an IQ-TREE Newick tree on a single outgroup while preserving
# branch lengths and support values, then ladderize clades in decreasing size.
#
# Also generates:
#   1. Rooted tree in Newick format
#   2. Equivalent rooted tree in NEXUS format
#   3. species_taxonomy.tsv for downstream orthology analyses
#
# Run from the repository root:
#   Rscript 04_orthology/07_root_tree.R

# ------------------------------- settings ------------------------------------
PROJECT  <- Sys.getenv("PROJECT", "results/04_orthology")
infile   <- file.path(PROJECT, "concat", "partitions_aa.txt.treefile")
outfile  <- file.path(PROJECT, "concat", "partitions_aa.txt_rooted_v2.treefile")
outgroup <- "Tribolium_castaneum"      # the tree is rooted on this tip
# -----------------------------------------------------------------------------

# 1) ape supplies the rooting and tree I/O.
suppressPackageStartupMessages(library(ape))

if (!file.exists(infile)) {
  stop("no such tree: ", infile, call. = FALSE)
}

# 2) Derive the NEXUS and taxonomy output paths from the tree path.
# Generate NEXUS filename from the requested Newick output filename
nexus_outfile <- sub(
  "\\.(treefile|nwk|tree|tre)$",
  ".nex",
  outfile,
  ignore.case = TRUE
)

# If the extension was not recognized, simply append ".nex"
if (identical(nexus_outfile, outfile)) {
  nexus_outfile <- paste0(outfile, ".nex")
}

# Assumes the rooted tree is being written inside PROJECT/concat/
# species_taxonomy.tsv will therefore be written to PROJECT/
project_dir <- dirname(dirname(outfile))

taxonomy_outfile <- file.path(
  project_dir,
  "species_taxonomy.tsv"
)

# 3) Read the unrooted ML tree from IQ-TREE.
tr <- read.tree(infile)

# 4) Check the outgroup is present as a tip.
# Require an exact tip-label match
if (!outgroup %in% tr$tip.label) {

  candidates <- grep(
    outgroup,
    tr$tip.label,
    value = TRUE,
    fixed = TRUE
  )

  msg <- paste0(
    "Outgroup '",
    outgroup,
    "' not found as an exact tip label."
  )

  if (length(candidates) > 0) {

    msg <- paste0(
      msg,
      "\nPossible matching label(s): ",
      paste(candidates, collapse = ", ")
    )
  }

  stop(msg, call. = FALSE)
}

# 5) Root on the outgroup, keeping branch lengths and support values.
# edgelabel = TRUE is important for IQ-TREE bootstrap/support labels because
# support values belong to branches and must remain attached to the same edges.

tr <- root(
  tr,
  outgroup = outgroup,
  resolve.root = TRUE,
  edgelabel = TRUE
)

# 6) Ladderize so clades are drawn in decreasing size.
# Arrange descendant clades in decreasing order
tr <- ladderize(
  tr,
  right = TRUE
)

# 7) Write the rooted tree in Newick.
# Branch lengths and support values are retained.
# The tree is NOT made ultrametric.

write.tree(
  tr,
  file = outfile
)

# 8) Write the same tree in NEXUS, which step 08 reads.
# Generate the Newick representation directly from the same R tree object.
# This ensures the NEXUS and .treefile versions contain the same topology,
# branch lengths, rooting, and support labels.

newick_string <- write.tree(tr)

cat(
  "#NEXUS\n",
  "begin trees;\n",
  "\ttree tree_1 = [&R] ",
  newick_string,
  "\n",
  "end;\n",
  sep = "",
  file = nexus_outfile
)

# 9) Map each tip to its family or superfamily.
# Explicit taxonomy for the species in the current OrthoFinder dataset.
#
# Species names must exactly match the names used in the phylogenetic tree
# and in the normalized OrthoFinder input.

taxonomy_map <- data.frame(

  Species = c(
    "Tribolium_castaneum",
    "Cylas_formicarius",
    "Sitophilus_oryzae",
    "Anthonomus_grandis",
    "Ceutorhynchus_assimilis",
    "Dendroctonus_ponderosae",
    "Hypothenemus_hampei",
    "Euwallacea_fornicatus",
    "Euwallacea_similis",
    "Exocentrus_adspersus",
    "Anoplophora_glabripennis",
    "Oncideres_impluviata",
    "Rhamnusium_bicolor",
    "Aromia_moschata",
    "Molorchus_minor",
    "Acanthoscelides_obtectus",
    "Phaedon_cochleariae",
    "Gonioctena_quinquepunctata",
    "Leptinotarsa_decemlineata",
    "Phyllotreta_striolata",
    "Psylliodes_chrysocephalus",
    "Diabrotica_balteata",
    "Diabrotica_virgifera",
    "Diorhabda_carinulata",
    "Diorhabda_sublineata"
  ),

  Family = c(
    "Tenebrionidae",
    "Brentidae",
    "Curculionidae",
    "Curculionidae",
    "Curculionidae",
    "Curculionidae",
    "Curculionidae",
    "Curculionidae",
    "Curculionidae",
    "Cerambycidae",
    "Cerambycidae",
    "Cerambycidae",
    "Cerambycidae",
    "Cerambycidae",
    "Cerambycidae",
    "Chrysomelidae",
    "Chrysomelidae",
    "Chrysomelidae",
    "Chrysomelidae",
    "Chrysomelidae",
    "Chrysomelidae",
    "Chrysomelidae",
    "Chrysomelidae",
    "Chrysomelidae",
    "Chrysomelidae"
  ),

  Superfamily = c(
    "Tenebrionoidea",
    "Curculionoidea",
    "Curculionoidea",
    "Curculionoidea",
    "Curculionoidea",
    "Curculionoidea",
    "Curculionoidea",
    "Curculionoidea",
    "Curculionoidea",
    "Chrysomeloidea",
    "Chrysomeloidea",
    "Chrysomeloidea",
    "Chrysomeloidea",
    "Chrysomeloidea",
    "Chrysomeloidea",
    "Chrysomeloidea",
    "Chrysomeloidea",
    "Chrysomeloidea",
    "Chrysomeloidea",
    "Chrysomeloidea",
    "Chrysomeloidea",
    "Chrysomeloidea",
    "Chrysomeloidea",
    "Chrysomeloidea",
    "Chrysomeloidea"
  ),

  stringsAsFactors = FALSE
)

# 10) Check every tip received an assignment.
# Every species in the tree must have a taxonomy entry.
missing_taxonomy <- setdiff(
  tr$tip.label,
  taxonomy_map$Species
)

if (length(missing_taxonomy) > 0) {

  stop(
    paste0(
      "The following tree species are missing from the taxonomy mapping:\n  ",
      paste(missing_taxonomy, collapse = "\n  ")
    ),
    call. = FALSE
  )
}

# Keep only species actually present in the tree.
# This also prevents stale taxonomy entries from entering downstream analyses.

taxonomy_tree <- taxonomy_map[
  taxonomy_map$Species %in% tr$tip.label,
]

# 11) Write species_taxonomy.tsv for the later steps.
write.table(
  taxonomy_tree,
  file = taxonomy_outfile,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# 12) Report what was written.
cat("\n")
cat("Input tree       :", infile, "\n")
cat("Outgroup         :", outgroup, "\n")
cat("Rooted           :", is.rooted(tr), "\n")
cat("Ultrametric      :", is.ultrametric(tr), "\n")
cat("Tips             :", Ntip(tr), "\n")
cat("Internal nodes   :", Nnode(tr), "\n")

cat(
  "Node labels      :",
  if (is.null(tr$node.label)) {
    0
  } else {
    sum(nzchar(tr$node.label))
  },
  "\n"
)

cat("\nOutputs:\n")
cat("Newick tree      :", outfile, "\n")
cat("NEXUS tree       :", nexus_outfile, "\n")
cat("Taxonomy table   :", taxonomy_outfile, "\n")
