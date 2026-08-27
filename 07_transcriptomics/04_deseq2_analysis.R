#!/usr/bin/env Rscript
#
# Differential expression with DESeq2.
#
# Salmon transcript abundances are summarised to gene level with tximport,
# size-factor normalised, and tested with the Wald test and Benjamini-Hochberg
# correction, at alpha = 0.05. The |log2FC| > 1 effect-size threshold used in
# the paper is applied downstream, in steps 05 and 06.
#
# Seven contrasts are run. Four have two or more biological replicates on each
# side and are reported as statistical results; three involve the larval
# samples, which come from a single individual, and are descriptive only. Each
# output file records which kind it is.
#
# Run from the repository root:
#   Rscript 07_transcriptomics/04_deseq2_analysis.R

# ------------------------------- settings ------------------------------------
SAMPLES_FILE <- Sys.getenv("SAMPLES_FILE", "07_transcriptomics/samples.tsv")
TX2GENE_FILE <- Sys.getenv("TX2GENE_FILE",
                           "results/07_transcriptomics/ref/tx2gene.tsv")
QUANT_DIR    <- Sys.getenv("QUANT_DIR", "results/07_transcriptomics/quant")
RESULTS_DIR  <- Sys.getenv("RESULTS_DIR", "results/07_transcriptomics/deseq2")
# -----------------------------------------------------------------------------

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# 1) Helper function to run DESeq2 comparison.
run_deseq_comparison <- function(dds_full, sample_ids, design_var, contrast, 
                                  comparison_name, results_dir, descriptive = FALSE) {
  #' Run a DESeq2 differential expression comparison
  #' 
  #' @param dds_full Full DESeqDataSet
  #' @param sample_ids Character vector of sample IDs to include
  #' @param design_var Name of the design variable (column in colData)
  #' @param contrast Vector for results() contrast argument
  #' @param comparison_name Name for output files
  #' @param results_dir Output directory
  #' @param descriptive If TRUE, adds warning to output about limited replicates
  
  message("\n", paste(rep("=", 70), collapse = ""))
  message("Comparison: ", comparison_name)
  if (descriptive) {
    message("*** DESCRIPTIVE ONLY - Limited biological replicates ***")
  }
  message(paste(rep("=", 70), collapse = ""))
  
  # Subset to relevant samples
  dds_sub <- dds_full[, colnames(dds_full) %in% sample_ids]
  
  # Show sample breakdown
  message("\nSample breakdown:")
  print(table(colData(dds_sub)[[design_var]]))
  
  # Drop unused factor levels and set design
  colData(dds_sub)[[design_var]] <- droplevels(colData(dds_sub)[[design_var]])
  design(dds_sub) <- as.formula(paste("~", design_var))
  
  # Run DESeq2
  message("\nRunning DESeq2...")
  dds_sub <- DESeq(dds_sub)
  
  # Extract results
  res <- results(dds_sub, contrast = contrast, alpha = 0.05)
  res <- res[order(res$padj), ]
  
  # Summary statistics
  n_sig_005 <- sum(res$padj < 0.05, na.rm = TRUE)
  n_sig_001 <- sum(res$padj < 0.01, na.rm = TRUE)
  n_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
  
  message("\nResults summary:")
  message("  Total genes tested: ", nrow(res))
  message("  Significant (padj < 0.05): ", n_sig_005)
  message("  Significant (padj < 0.01): ", n_sig_001)
  message("  Up-regulated (", contrast[2], " > ", contrast[3], "): ", n_up)
  message("  Down-regulated (", contrast[2], " < ", contrast[3], "): ", n_down)
  
  # Create output dataframe with metadata
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE", 
                        "stat", "pvalue", "padj")]
  
  # Add comparison metadata
  attr(res_df, "comparison") <- comparison_name
  attr(res_df, "contrast") <- paste(contrast, collapse = " vs ")
  attr(res_df, "descriptive") <- descriptive
  attr(res_df, "n_samples") <- paste(table(colData(dds_sub)[[design_var]]), collapse = " vs ")
  
  # Save results
  outfile <- file.path(results_dir, paste0("DE_", comparison_name, ".csv"))
  
  # Add header comment for descriptive analyses
  if (descriptive) {
    header <- paste0(
      "# DESCRIPTIVE ANALYSIS - LIMITED BIOLOGICAL REPLICATES\n",
      "# Comparison: ", comparison_name, "\n",
      "# Contrast: ", paste(contrast, collapse = " vs "), "\n",
      "# WARNING: P-values should be interpreted with caution.\n",
      "# Larval samples (L1He, L1Mg) are from a single individual.\n",
      "# Focus on large fold-changes for biological interpretation.\n"
    )
    writeLines(header, outfile)
    write.table(res_df, outfile, sep = ",", row.names = FALSE, 
                col.names = TRUE, append = TRUE)
  } else {
    write.csv(res_df, outfile, row.names = FALSE)
  }
  
  message("  Saved: ", basename(outfile))
  
  return(list(dds = dds_sub, results = res, results_df = res_df))
}

# 2) Load packages.
message("\n", paste(rep("=", 70), collapse = ""))
message("DESeq2 Analysis Pipeline - Oncideres impluviata")
message(paste(rep("=", 70), collapse = ""))
message("\nLoading packages...")

suppressPackageStartupMessages({
  library(tximport)
  library(DESeq2)
  library(tidyverse)
})

message("  DESeq2 version: ", packageVersion("DESeq2"))
message("  tximport version: ", packageVersion("tximport"))

# 3) Load sample metadata.
message("\n[Step 1] Loading sample metadata...")

samples <- read_tsv(SAMPLES_FILE, show_col_types = FALSE)
message("  Total samples: ", nrow(samples))

# Create file paths
samples$files <- file.path(QUANT_DIR, samples$sample_id, "quant.sf")

# Verify all quant files exist
missing_files <- samples$files[!file.exists(samples$files)]
if (length(missing_files) > 0) {
  stop("Missing quant.sf files:\n  ", paste(missing_files, collapse = "\n  "))
}
message("  All quant.sf files found")

# Display sample summary
message("\nSample breakdown:")
message("  By stage:")
print(table(samples$stage))
message("  By stage and sex:")
print(table(samples$stage, samples$sex))
message("  By stage, sex, and structure:")
print(table(samples$stage, samples$sex, samples$structure))

# Save sample summary to file
summary_file <- file.path(RESULTS_DIR, "sample_summary.txt")
sink(summary_file)
cat("Sample Summary for DESeq2 Analysis\n")
cat("===================================\n\n")
cat("Total samples:", nrow(samples), "\n\n")
cat("By stage:\n")
print(table(samples$stage))
cat("\nBy stage and sex:\n")
print(table(samples$stage, samples$sex))
cat("\nDetailed sample list:\n")
print(samples %>% select(sample_id, stage, sex, structure))
cat("\n\nComparison Design:\n")
cat("==================\n")
cat("\nSTATISTICALLY VALID:\n")
cat("1. Male vs Female (antennae): AM1An, AM2An vs AF1An, AF2An (2 vs 2)\n")
cat("2. Male vs Female (thorax): AM2Th vs AF1Th, AF2Th (1 vs 2) - LIMITED\n")
cat("3. Male vs Female (pooled): 3 male vs 4 female adults\n")
cat("4. Antennae vs Thorax: 4 antennae vs 3 thorax\n")
cat("\nDESCRIPTIVE ONLY (single larval individual):\n")
cat("5. Larva vs Adult Male: L1He, L1Mg vs 3 male samples\n")
cat("6. Larva vs Adult Female: L1He, L1Mg vs 4 female samples\n")
cat("7. Larva vs All Adults: L1He, L1Mg vs 7 adult samples\n")
sink()
message("  Saved: ", basename(summary_file))

# 4) Load transcript-to-gene mapping.
message("\n[Step 2] Loading tx2gene mapping...")
tx2gene <- read_tsv(TX2GENE_FILE, show_col_types = FALSE)
message("  Transcripts: ", nrow(tx2gene))
message("  Unique genes: ", n_distinct(tx2gene$gene_id))

# 5) Import salmon counts with tximport.
message("\n[Step 3] Importing Salmon quantifications...")

txi <- tximport(
  samples$files, 
  type = "salmon", 
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE,
  countsFromAbundance = "lengthScaledTPM"
)

colnames(txi$counts) <- samples$sample_id
colnames(txi$abundance) <- samples$sample_id
colnames(txi$length) <- samples$sample_id

message("  Genes imported: ", nrow(txi$counts))
message("  Samples: ", ncol(txi$counts))

# 6) Create master DESeq2 dataset.
message("\n[Step 4] Creating DESeq2 dataset...")

# Prepare sample metadata
sample_info <- samples %>%
  column_to_rownames("sample_id") %>%
  dplyr::select(stage, sex, structure, group, batch) %>%
  mutate(
    stage = factor(stage, levels = c("larva", "adult")),
    sex = factor(sex),
    structure = factor(structure),
    # Create combined condition for flexible subsetting
    stage_sex = factor(paste(stage, sex, sep = "_")),
    condition = factor(paste(stage, sex, structure, sep = "_"))
  )

# Create DESeqDataSet with general design (will be updated per comparison)
dds <- DESeqDataSetFromTximport(
  txi, 
  colData = sample_info, 
  design = ~ condition
)

message("  Initial genes: ", nrow(dds))

# Filter lowly expressed genes
# Keep genes with >= 10 counts in at least 2 samples
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]
message("  Genes after filtering (>=10 counts in >=2 samples): ", nrow(dds))

# Run DESeq on full dataset (for normalized counts)
message("\n[Step 5] Running DESeq2 on full dataset...")
dds <- DESeq(dds)

# 7) Generate normalized counts (for visualization).
message("\n[Step 6] Generating normalized counts...")

# Variance stabilizing transformation (for heatmaps, PCA)
vsd <- vst(dds, blind = FALSE)
norm_counts_vst <- assay(vsd)

# DESeq2 size-factor normalized counts
norm_counts_sf <- counts(dds, normalized = TRUE)

# Save normalized counts
write.csv(norm_counts_vst, file.path(RESULTS_DIR, "normalized_counts_vst.csv"))
write.csv(norm_counts_sf, file.path(RESULTS_DIR, "normalized_counts_deseq2.csv"))
message("  Saved: normalized_counts_vst.csv")
message("  Saved: normalized_counts_deseq2.csv")

# 8) Define sample groups for each comparison.
# Adult samples
adult_male_ant <- samples %>% filter(stage == "adult", sex == "male", structure == "antenna") %>% pull(sample_id)
adult_female_ant <- samples %>% filter(stage == "adult", sex == "female", structure == "antenna") %>% pull(sample_id)
adult_male_th <- samples %>% filter(stage == "adult", sex == "male", structure == "thorax") %>% pull(sample_id)
adult_female_th <- samples %>% filter(stage == "adult", sex == "female", structure == "thorax") %>% pull(sample_id)

# Combined groups
all_adult_male <- c(adult_male_ant, adult_male_th)
all_adult_female <- c(adult_female_ant, adult_female_th)
all_adults <- c(all_adult_male, all_adult_female)
all_antennae <- c(adult_male_ant, adult_female_ant)
all_thorax <- c(adult_male_th, adult_female_th)

# Larval samples
larva_samples <- samples %>% filter(stage == "larva") %>% pull(sample_id)

message("\nSample groups defined:")
message("  Adult male antennae: ", paste(adult_male_ant, collapse = ", "))
message("  Adult female antennae: ", paste(adult_female_ant, collapse = ", "))
message("  Adult male thorax: ", paste(adult_male_th, collapse = ", "))
message("  Adult female thorax: ", paste(adult_female_th, collapse = ", "))
message("  Larva: ", paste(larva_samples, collapse = ", "))

# Store all results
all_results <- list()

# 9) Statistically valid comparisons.
message("\n", paste(rep("#", 70), collapse = ""))
message("# STATISTICALLY VALID COMPARISONS")
message(paste(rep("#", 70), collapse = ""))

# 10) 1: adult male vs female (antennae only) - 2 vs 2.
all_results$antennae_MvF <- run_deseq_comparison(
  dds_full = dds,
  sample_ids = c(adult_male_ant, adult_female_ant),
  design_var = "sex",
  contrast = c("sex", "male", "female"),
  comparison_name = "antennae_male_vs_female",
  results_dir = RESULTS_DIR,
  descriptive = FALSE
)

# 11) 2: adult male vs female (thorax) - 1 vs 2 (limited).
all_results$thorax_MvF <- run_deseq_comparison(
  dds_full = dds,
  sample_ids = c(adult_male_th, adult_female_th),
  design_var = "sex",
  contrast = c("sex", "male", "female"),
  comparison_name = "thorax_male_vs_female",
  results_dir = RESULTS_DIR,
  descriptive = FALSE  # Technically runs but limited power
)

# 12) 3: adult male vs female (all tissues pooled) - 3 vs 4 new.
all_results$adults_MvF <- run_deseq_comparison(
  dds_full = dds,
  sample_ids = all_adults,
  design_var = "sex",
  contrast = c("sex", "male", "female"),
  comparison_name = "adults_pooled_male_vs_female",
  results_dir = RESULTS_DIR,
  descriptive = FALSE
)

# 13) 4: antennae vs thorax (all adults) - 4 vs 3.
all_results$tissue <- run_deseq_comparison(
  dds_full = dds,
  sample_ids = all_adults,
  design_var = "structure",
  contrast = c("structure", "antenna", "thorax"),
  comparison_name = "antennae_vs_thorax",
  results_dir = RESULTS_DIR,
  descriptive = FALSE
)

# 14) Descriptive comparisons (limited biological replicates).
message("\n", paste(rep("#", 70), collapse = ""))
message("# DESCRIPTIVE COMPARISONS - Interpret with caution")
message("# Larval samples are from a single individual")
message(paste(rep("#", 70), collapse = ""))

# 15) 5: larva vs adult male - 2 vs 3 new (descriptive).
all_results$larva_vs_male <- run_deseq_comparison(
  dds_full = dds,
  sample_ids = c(larva_samples, all_adult_male),
  design_var = "stage",
  contrast = c("stage", "larva", "adult"),
  comparison_name = "larva_vs_adult_male",
  results_dir = RESULTS_DIR,
  descriptive = TRUE
)

# 16) 6: larva vs adult female - 2 vs 4 new (descriptive).
all_results$larva_vs_female <- run_deseq_comparison(
  dds_full = dds,
  sample_ids = c(larva_samples, all_adult_female),
  design_var = "stage",
  contrast = c("stage", "larva", "adult"),
  comparison_name = "larva_vs_adult_female",
  results_dir = RESULTS_DIR,
  descriptive = TRUE
)

# 17) 7: larva vs all adults - 2 vs 7 new (descriptive).
all_results$larva_vs_adults <- run_deseq_comparison(
  dds_full = dds,
  sample_ids = c(larva_samples, all_adults),
  design_var = "stage",
  contrast = c("stage", "larva", "adult"),
  comparison_name = "larva_vs_all_adults",
  results_dir = RESULTS_DIR,
  descriptive = TRUE
)

# 18) Save all R objects.
message("\n", paste(rep("=", 70), collapse = ""))
message("Saving R objects...")
message(paste(rep("=", 70), collapse = ""))

save(
  dds,              # Full DESeqDataSet
  vsd,              # VST-transformed data
  norm_counts_vst,  # VST normalized counts matrix
  norm_counts_sf,   # Size-factor normalized counts
  all_results,      # All comparison results
  samples,          # Sample metadata
  sample_info,      # DESeq2 colData
  tx2gene,          # Transcript to gene mapping
  file = file.path(RESULTS_DIR, "deseq2_results.RData")
)

message("  Saved: deseq2_results.RData")

# 19) Final summary.
message("\n", paste(rep("=", 70), collapse = ""))
message("ANALYSIS COMPLETE")
message(paste(rep("=", 70), collapse = ""))

message("\nResults directory: ", RESULTS_DIR)
message("\nOutput files:")
message("  - normalized_counts_vst.csv")
message("  - normalized_counts_deseq2.csv")
message("  - DE_antennae_male_vs_female.csv")
message("  - DE_thorax_male_vs_female.csv")
message("  - DE_adults_pooled_male_vs_female.csv")
message("  - DE_antennae_vs_thorax.csv")
message("  - DE_larva_vs_adult_male.csv (DESCRIPTIVE)")
message("  - DE_larva_vs_adult_female.csv (DESCRIPTIVE)")
message("  - DE_larva_vs_all_adults.csv (DESCRIPTIVE)")
message("  - deseq2_results.RData")
message("  - sample_summary.txt")

message("\n", paste(rep("-", 70), collapse = ""))
message("Summary of significant genes (padj < 0.05):")
message(paste(rep("-", 70), collapse = ""))

for (name in names(all_results)) {
  res <- all_results[[name]]$results
  n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
  is_desc <- ifelse(grepl("larva", name), " [DESCRIPTIVE]", "")
  message(sprintf("  %-30s: %5d genes%s", name, n_sig, is_desc))
}

message("\n", paste(rep("=", 70), collapse = ""))
message("Pipeline completed: ", Sys.time())
message(paste(rep("=", 70), collapse = ""))
