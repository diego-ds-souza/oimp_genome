#!/usr/bin/env Rscript
#
# Are differentially expressed genes overrepresented inside tandem clusters?
#
# Every clustered gene is by construction a member of one of the 16
# host-interface families, and those families are differentially expressed more
# often than the average annotated gene. A genome-wide 2x2 table therefore
# compares two groups that differ in family membership as well as in
# clustering, which is what the reviewer questioned.
#
# The script runs the genome-wide test, then repeats it restricted to the 16
# target families - so clustered genes are compared with unclustered members of
# the same families - and separately for stage-biased and for tissue- or
# sex-biased genes. It also reports a Mantel-Haenszel estimate stratified by
# family, which controls for family size and identity.
#
# The odds ratio reported is the sample odds ratio (a*d)/(b*c), not the
# conditional maximum-likelihood estimate that fisher.test() returns.
#
# Run from the repository root:
#   Rscript 08_chromosomal_organization/03_cluster_de_enrichment.R

# ------------------------------- settings ------------------------------------
DE_DIR      <- Sys.getenv("DE_DIR", "results/07_transcriptomics/deseq2")
CHR_DIR     <- Sys.getenv("CHR_DIR", "results/08_chromosomal_organization")
OUT_DIR     <- Sys.getenv("OUT_DIR", "results/08_chromosomal_organization/enrichment")
POSITIONS   <- file.path(CHR_DIR, "gene_positions_chromosome_level.tsv")
CLUSTERS    <- file.path(CHR_DIR, "gene_clusters_chromosome_level.tsv")
PADJ_CUTOFF <- 0.05      # adjusted p-value threshold
LFC_CUTOFF  <- 1         # |log2 fold change| threshold
STAGE_PREFIX <- "larva"  # DE_*.csv files starting with this are stage contrasts
# -----------------------------------------------------------------------------

# Contrasts whose DE_*.csv file name begins with this prefix are developmental
# (stage) contrasts; all others are adult tissue or sex contrasts.
STAGE_PREFIX <- "larva"

# When a gene carries two family assignments it is collapsed to a single primary
# family before stratification, using the same precedence as Tables S16-S18
# (UGT and SerPro take precedence over CAZyme; OR over GR). Lower rank wins.
FAMILY_RANK <- c(CAZyme = 2, GR = 1)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
LOG   <- file.path(OUT_DIR, paste0("cluster_DE_enrichment_", STAMP, ".log"))

log_msg <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n", sep = "")
  cat(msg, "\n", sep = "", file = LOG, append = TRUE)
}

log_msg(strrep("=", 78))
log_msg("Tandem cluster enrichment of differentially expressed genes")
log_msg(strrep("=", 78))
log_msg("DE results   : ", DE_DIR)
log_msg("Positions    : ", POSITIONS)
log_msg("Clusters     : ", CLUSTERS)
log_msg("Thresholds   : padj < ", PADJ_CUTOFF, " and |log2FC| > ", LFC_CUTOFF)
log_msg("Started      : ", format(Sys.time()))
log_msg("")

for (f in c(POSITIONS, CLUSTERS)) {
  if (!file.exists(f)) stop("Input not found: ", f)
}

# 1) Gene universe, family assignment and cluster membership.
log_msg("[Step 1] Reading gene positions and cluster membership...")

positions <- read.delim(POSITIONS, stringsAsFactors = FALSE)
universe  <- unique(positions$gene_id)

assigned <- positions[!positions$family %in% c("NA", "", NA), c("gene_id", "family")]
assigned$rank <- ifelse(assigned$family %in% names(FAMILY_RANK),
                        FAMILY_RANK[assigned$family], 0)
assigned <- assigned[order(assigned$gene_id, assigned$rank), ]
primary  <- assigned[!duplicated(assigned$gene_id), ]

family  <- setNames(primary$family, primary$gene_id)
targets <- primary$gene_id

clusters  <- read.delim(CLUSTERS, stringsAsFactors = FALSE)
clustered <- unique(unlist(strsplit(clusters$genes, ",", fixed = TRUE)))

log_msg("  Annotated genes with chromosome-level positions : ", length(universe))
log_msg("  Genes assigned to a target family               : ", length(targets))
log_msg("  Tandem clusters                                 : ", nrow(clusters))
log_msg("  Genes in tandem clusters                        : ", length(clustered))
log_msg("")

# 2) Differential expression calls, per contrast.
log_msg("[Step 2] Reading DESeq2 results...")

de_files <- sort(list.files(DE_DIR, pattern = "^DE_.*\\.csv$", full.names = TRUE))
if (length(de_files) == 0) stop("No DE_*.csv files found in ", DE_DIR)

de_sets <- list()
for (f in de_files) {
  contrast <- sub("^DE_", "", file_path_sans_ext(basename(f)))
  d <- read.csv(f, stringsAsFactors = FALSE, comment.char = "#")
  keep <- !is.na(d$padj) & !is.na(d$log2FoldChange) &
          d$padj < PADJ_CUTOFF & abs(d$log2FoldChange) > LFC_CUTOFF
  de_sets[[contrast]] <- unique(d$gene_id[keep])
  log_msg(sprintf("  %-34s %6d tested  %6d DE", contrast, nrow(d), sum(keep)))
}

is_stage    <- startsWith(names(de_sets), STAGE_PREFIX)
de_any      <- unique(unlist(de_sets))
de_stage    <- unique(unlist(de_sets[is_stage]))
de_tissuesex<- unique(unlist(de_sets[!is_stage]))

log_msg("")
log_msg("  Stage contrasts      : ", paste(names(de_sets)[is_stage], collapse = ", "))
log_msg("  Tissue/sex contrasts : ", paste(names(de_sets)[!is_stage], collapse = ", "))
log_msg("")

# 3) Tests.
two_by_two <- function(u, de) {
  in_cl <- u %in% clustered
  in_de <- u %in% de
  c(a = sum( in_cl &  in_de),   # clustered, differentially expressed
    b = sum( in_cl & !in_de),   # clustered, not
    c = sum(!in_cl &  in_de),   # unclustered, differentially expressed
    d = sum(!in_cl & !in_de))   # unclustered, not
}

run_test <- function(label, u, de, stratify) {
  n <- two_by_two(u, de)
  a <- n[["a"]]; b <- n[["b"]]; cc <- n[["c"]]; d <- n[["d"]]
  or <- (a * d) / (b * cc)                       # sample odds ratio
  p  <- fisher.test(matrix(c(a, b, cc, d), nrow = 2, byrow = TRUE))$p.value
  mh <- NA_real_
  if (stratify) {
    num <- 0; den <- 0
    for (fm in unique(family)) {
      m  <- u[family[u] %in% fm]
      nn <- two_by_two(m, de)
      tt <- sum(nn)
      if (tt > 0) {
        num <- num + nn[["a"]] * nn[["d"]] / tt
        den <- den + nn[["b"]] * nn[["c"]] / tt
      }
    }
    mh <- if (den > 0) num / den else NA_real_
  }
  data.frame(test = label, a = a, b = b, c = cc, d = d,
             odds_ratio = or, p_value = p, mantel_haenszel_or = mh,
             stringsAsFactors = FALSE)
}

log_msg("[Step 3] Running tests...")
log_msg("")

results <- rbind(
  run_test("genome-wide, DE in any contrast",
           universe, de_any, stratify = FALSE),
  run_test("within target families, DE in any contrast",
           targets, de_any, stratify = TRUE),
  run_test("within target families, tissue- or sex-biased",
           targets, de_tissuesex, stratify = TRUE),
  run_test("within target families, stage-biased",
           targets, de_stage, stratify = TRUE)
)

log_msg(sprintf("  %-46s %6s %6s %7s %7s %8s %11s %8s",
                "test", "a", "b", "c", "d", "OR", "p", "MH OR"))
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  log_msg(sprintf("  %-46s %6d %6d %7d %7d %8.2f %11.2e %8s",
                  r$test, r$a, r$b, r$c, r$d, r$odds_ratio, r$p_value,
                  ifelse(is.na(r$mantel_haenszel_or), "-",
                         sprintf("%.2f", r$mantel_haenszel_or))))
}

out_tsv <- file.path(OUT_DIR, "cluster_DE_enrichment.tsv")
write.table(results, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

log_msg("")
log_msg("  a = clustered and differentially expressed")
log_msg("  b = clustered, not differentially expressed")
log_msg("  c = unclustered and differentially expressed")
log_msg("  d = unclustered, not differentially expressed")
log_msg("  Odds ratio is (a*d)/(b*c); p from two-sided fisher.test().")
log_msg("  Mantel-Haenszel odds ratio is pooled over strata defined by primary gene family.")
log_msg("")
log_msg("Written: ", out_tsv)
log_msg("Written: ", LOG)
log_msg("Completed: ", format(Sys.time()))
