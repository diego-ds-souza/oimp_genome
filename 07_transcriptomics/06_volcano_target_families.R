#!/usr/bin/env Rscript
#
# Figure 6: volcano plots of the host-interface gene families.
#
# DESeq2 results are filtered to genes assigned to a chemosensory,
# detoxification or digestion family, then plotted per functional group, as a
# combined panel, and faceted for the antennae-versus-thorax contrast. The top
# genes per family are labelled and written out for identification.
#
# The same thresholds as the genome-wide analysis are applied here: adjusted
# p < 0.05 and |log2FC| > 1.
#
# Run from the repository root:
#   Rscript 07_transcriptomics/06_volcano_target_families.R

# ------------------------------- settings ------------------------------------
DESEQ_RESULTS   <- Sys.getenv("DESEQ_RESULTS",
                              "results/07_transcriptomics/deseq2/deseq2_results.RData")
FAMILIES_DIR    <- Sys.getenv("FAMILIES_DIR", "results/07_transcriptomics/families")
GENE_NAMES_FILE <- Sys.getenv("GENE_NAMES_FILE",
                              "results/07_transcriptomics/ref/gene_names.tsv")
OUTPUT_DIR      <- Sys.getenv("OUTPUT_DIR",
                              "results/07_transcriptomics/volcano_target_families")
PLOTS_DIR       <- file.path(OUTPUT_DIR, "plots")
PADJ_CUTOFF     <- 0.05     # adjusted p-value threshold
LFC_CUTOFF      <- 1        # |log2 fold change| threshold
TOP_N_LABEL     <- 10       # top genes labelled per family
# -----------------------------------------------------------------------------

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

# 1) Load packages.
message("\n", paste(rep("=", 70), collapse = ""))
message("Family-Specific Volcano Plots")
message(paste(rep("=", 70), collapse = ""))
message("\nLoading packages...")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) library(patchwork)

# 2) Load data.
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 1: Loading Data")
message(paste(rep("=", 70), collapse = ""))

# Load DESeq2 results
message("\n[Step 1] Loading DESeq2 results...")
load(DESEQ_RESULTS)
message("  Comparisons: ", paste(names(all_results), collapse = ", "))

# Load gene names
message("\n[Step 2] Loading gene annotations...")
gene_names <- read_tsv(GENE_NAMES_FILE, show_col_types = FALSE)
message("  Genes: ", nrow(gene_names))

# 3) Load gene family assignments.
message("\n[Step 3] Loading gene family assignments...")

FAMILY_GROUPS <- list(
  Chemosensory = c("OBP", "CSP", "OR", "GR", "IR", "SNMP"),
  Detoxification = c("P450", "GST", "UGT", "CCE", "ABC"),
  Digestion = c("SerPro", "CysPro", "AspPro", "MetPro", "CAZyme")
)

get_family_genes <- function(family_name, families_dir) {
  fam_file <- file.path(families_dir, paste0("genes_", family_name, ".txt"))
  if (file.exists(fam_file)) {
    genes <- readLines(fam_file)
    return(genes[genes != ""])
  }
  return(character(0))
}

family_assignments <- data.frame()
for (group_name in names(FAMILY_GROUPS)) {
  for (fam_name in FAMILY_GROUPS[[group_name]]) {
    genes <- get_family_genes(fam_name, FAMILIES_DIR)
    if (length(genes) > 0) {
      family_assignments <- rbind(family_assignments, data.frame(
        gene_id = genes,
        family = fam_name,
        group = group_name,
        stringsAsFactors = FALSE
      ))
      message("  ", fam_name, ": ", length(genes), " genes")
    }
  }
}

message("  Total genes in target families: ", nrow(family_assignments))

# Add gene names
family_assignments <- family_assignments %>%
  left_join(gene_names %>% select(gene_id, name, product), by = "gene_id")

# 4) Define colors and labels.
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 2: Setting Up Plot Parameters")
message(paste(rep("=", 70), collapse = ""))

# Family colors
family_colors <- c(
  # Chemosensory - greens
  "OBP" = "#1B9E77",
  "CSP" = "#66C2A5", 
  "OR" = "#D95F02",   # Orange for visibility
  "GR" = "#E7298A",
  "IR" = "#7570B3",
  "SNMP" = "#A6761D",
  # Detoxification - oranges/reds
  "P450" = "#E6AB02",
  "GST" = "#666666",
  "UGT" = "#1F78B4",
  "CCE" = "#E31A1C",
  "ABC" = "#FB9A99",
  # Digestion - blues
  "SerPro" = "#33A02C",
  "CysPro" = "#6A3D9A",
  "AspPro" = "#A6CEE3",
  "MetPro" = "#CAB2D6",
  "CAZyme" = "#B15928"
)

# Group colors (for group-level plots)
group_colors <- c(
  "Chemosensory" = "#1B9E77",
  "Detoxification" = "#D95F02",
  "Digestion" = "#1F78B4"
)

# Comparison labels
comparison_labels <- c(
  "antennae_MvF" = "Male vs Female\n(Antennae)",
  "thorax_MvF" = "Male vs Female\n(Thorax)",
  "adults_MvF" = "Male vs Female\n(All Adults)",
  "tissue" = "Antennae vs Thorax\n(Adults)",
  "larva_vs_male" = "Larva vs Male",
  "larva_vs_female" = "Larva vs Female",
  "larva_vs_adults" = "Larva vs Adults\n(All)"
)

# Key comparisons for main figures
key_comparisons <- c("larva_vs_adults", "tissue")

# 5) Prepare DE data for target families.
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 3: Preparing DE Data for Target Families")
message(paste(rep("=", 70), collapse = ""))

message("\n[Step 4] Extracting DE results for target family genes...")

all_de_data <- data.frame()

for (comp_name in names(all_results)) {
  
  # Get DE results
  de_res <- as.data.frame(all_results[[comp_name]]$results)
  de_res$gene_id <- rownames(de_res)
  
  # Filter to target family genes only
  de_family <- de_res %>%
    filter(gene_id %in% family_assignments$gene_id) %>%
    left_join(family_assignments, by = "gene_id") %>%
    mutate(
      comparison = comp_name,
      comparison_label = comparison_labels[comp_name],
      # Significance categories
      sig_status = case_when(
        is.na(padj) ~ "NS",
        padj >= PADJ_CUTOFF ~ "NS",
        abs(log2FoldChange) < LFC_CUTOFF ~ "NS",
        log2FoldChange > 0 ~ "Up",
        log2FoldChange < 0 ~ "Down"
      ),
      # For plotting
      neg_log10_padj = -log10(pmax(padj, 1e-300))  # Avoid -Inf
    )
  
  all_de_data <- rbind(all_de_data, de_family)
}

message("  Total data points: ", nrow(all_de_data))

# Summary
de_summary <- all_de_data %>%
  filter(sig_status != "NS") %>%
  group_by(comparison, group, family) %>%
  summarise(
    n_up = sum(sig_status == "Up"),
    n_down = sum(sig_status == "Down"),
    n_total = n(),
    .groups = "drop"
  )

message("\n  DE genes per family (Larva vs Adults):")
print(
  de_summary %>% 
    filter(comparison == "larva_vs_adults") %>%
    arrange(group, desc(n_total)) %>%
    select(group, family, n_up, n_down, n_total)
)

# 6) Create volcano plots by functional group.
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 4: Creating Volcano Plots by Functional Group")
message(paste(rep("=", 70), collapse = ""))

message("\n[Step 5] Creating group-specific volcano plots...")

for (comp_name in key_comparisons) {
  
  comp_label <- comparison_labels[comp_name]
  comp_data <- all_de_data %>% filter(comparison == comp_name)
  
  for (group_name in names(FAMILY_GROUPS)) {
    
    group_data <- comp_data %>% filter(group == group_name)
    
    if (nrow(group_data) == 0) next
    
    # Identify top genes to label (top N per family by padj)
    top_genes <- group_data %>%
      filter(sig_status != "NS") %>%
      group_by(family) %>%
      slice_min(padj, n = TOP_N_LABEL, with_ties = FALSE) %>%
      ungroup()
    
    # Create label column
    group_data <- group_data %>%
      mutate(
        label = ifelse(gene_id %in% top_genes$gene_id, name, NA_character_)
      )
    
    # Calculate axis limits
    max_lfc <- max(abs(group_data$log2FoldChange), na.rm = TRUE) * 1.1
    max_pval <- max(group_data$neg_log10_padj[is.finite(group_data$neg_log10_padj)], na.rm = TRUE) * 1.1
    
    # Create plot
    p <- ggplot(group_data, aes(x = log2FoldChange, y = neg_log10_padj)) +
      # Non-significant points (grey)
      geom_point(
        data = group_data %>% filter(sig_status == "NS"),
        color = "grey70", alpha = 0.5, size = 1.5
      ) +
      # Significant points (colored by family)
      geom_point(
        data = group_data %>% filter(sig_status != "NS"),
        aes(color = family),
        alpha = 0.8, size = 2.5
      ) +
      # Labels for top genes
      geom_text_repel(
        data = group_data %>% filter(!is.na(label)),
        aes(label = label, color = family),
        size = 3,
        max.overlaps = 20,
        segment.color = "grey50",
        segment.alpha = 0.5,
        fontface = "bold",
        show.legend = FALSE
      ) +
      # Significance thresholds
      geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", color = "grey40") +
      # Scales
      scale_color_manual(values = family_colors, name = "Gene Family") +
      scale_x_continuous(limits = c(-max_lfc, max_lfc)) +
      coord_cartesian(ylim = c(0, max_pval)) +
      # Labels
      labs(
        title = paste0(group_name, " Genes — ", gsub("\n", " ", comp_label)),
        subtitle = paste0(
          sum(group_data$sig_status == "Up"), " up, ",
          sum(group_data$sig_status == "Down"), " down (padj < ", PADJ_CUTOFF, ", |log2FC| > ", LFC_CUTOFF, ")"
        ),
        x = expression(log[2]~"Fold Change"),
        y = expression(-log[10]~"adjusted p-value")
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "grey40"),
        legend.position = "right",
        panel.grid.minor = element_blank()
      )
    
    # Save
    filename <- paste0("volcano_", group_name, "_", comp_name, ".pdf")
    ggsave(file.path(PLOTS_DIR, filename), p, width = 12, height = 8)
    message("  Saved: ", filename)
  }
}

# 7) Create combined volcano plots (all families).
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 5: Creating Combined Volcano Plots")
message(paste(rep("=", 70), collapse = ""))

message("\n[Step 6] Creating combined volcano plots with all families...")

for (comp_name in key_comparisons) {
  
  comp_label <- comparison_labels[comp_name]
  comp_data <- all_de_data %>% filter(comparison == comp_name)
  
  # Identify top genes to label (top 3 per family)
  top_genes <- comp_data %>%
    filter(sig_status != "NS") %>%
    group_by(family) %>%
    slice_min(padj, n = 3, with_ties = FALSE) %>%
    ungroup()
  
  comp_data <- comp_data %>%
    mutate(
      label = ifelse(gene_id %in% top_genes$gene_id, name, NA_character_)
    )
  
  max_lfc <- max(abs(comp_data$log2FoldChange), na.rm = TRUE) * 1.1
  max_pval <- max(comp_data$neg_log10_padj[is.finite(comp_data$neg_log10_padj)], na.rm = TRUE) * 1.1
  
  p <- ggplot(comp_data, aes(x = log2FoldChange, y = neg_log10_padj)) +
    # Non-significant
    geom_point(
      data = comp_data %>% filter(sig_status == "NS"),
      color = "grey80", alpha = 0.4, size = 1
    ) +
    # Significant - colored by group
    geom_point(
      data = comp_data %>% filter(sig_status != "NS"),
      aes(color = group),
      alpha = 0.7, size = 2
    ) +
    # Labels
    geom_text_repel(
      data = comp_data %>% filter(!is.na(label)),
      aes(label = label, color = group),
      size = 2.5,
      max.overlaps = 30,
      segment.color = "grey50",
      segment.alpha = 0.3,
      show.legend = FALSE
    ) +
    # Thresholds
    geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", color = "grey40") +
    # Scales
    scale_color_manual(values = group_colors, name = "Functional Group") +
    scale_x_continuous(limits = c(-max_lfc, max_lfc)) +
    coord_cartesian(ylim = c(0, max_pval)) +
    labs(
      title = paste0("Target Gene Families — ", gsub("\n", " ", comp_label)),
      subtitle = paste0(
        "Chemosensory + Detoxification + Digestion genes | ",
        sum(comp_data$sig_status != "NS"), " significantly DE"
      ),
      x = expression(log[2]~"Fold Change"),
      y = expression(-log[10]~"adjusted p-value")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      legend.position = "right"
    )
  
  filename <- paste0("volcano_all_families_", comp_name, ".pdf")
  ggsave(file.path(PLOTS_DIR, filename), p, width = 12, height = 8)
  message("  Saved: ", filename)
}

# 8) Create faceted volcano plot (all groups, key comparison).
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 6: Creating Faceted Volcano Plot")
message(paste(rep("=", 70), collapse = ""))

message("\n[Step 7] Creating faceted volcano plot (Larva vs Adults)...")

comp_data <- all_de_data %>% filter(comparison == "larva_vs_adults")

# Top 5 genes per group
top_genes <- comp_data %>%
  filter(sig_status != "NS") %>%
  group_by(group) %>%
  slice_min(padj, n = 5, with_ties = FALSE) %>%
  ungroup()

comp_data <- comp_data %>%
  mutate(label = ifelse(gene_id %in% top_genes$gene_id, name, NA_character_))

p_facet <- ggplot(comp_data, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(
    data = comp_data %>% filter(sig_status == "NS"),
    color = "grey80", alpha = 0.4, size = 1
  ) +
  geom_point(
    data = comp_data %>% filter(sig_status != "NS"),
    aes(color = family),
    alpha = 0.8, size = 2
  ) +
  geom_text_repel(
    data = comp_data %>% filter(!is.na(label)),
    aes(label = label),
    size = 2.5,
    max.overlaps = 15,
    segment.color = "grey50"
  ) +
  geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = family_colors, name = "Gene Family") +
  facet_wrap(~ group, scales = "free") +
  labs(
    title = "Target Gene Families — Larva vs Adults",
    subtitle = "Faceted by functional group | Top 5 genes labeled per group",
    x = expression(log[2]~"Fold Change"),
    y = expression(-log[10]~"adjusted p-value")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2))

ggsave(file.path(PLOTS_DIR, "volcano_faceted_larva_vs_adults.pdf"), 
       p_facet, width = 14, height = 6)
message("  Saved: volcano_faceted_larva_vs_adults.pdf")

# 9) Faceted volcano plot for antennae vs thorax.
message("\n[Step 7b] Creating faceted volcano plot (Antennae vs Thorax)...")

comp_data_tissue <- all_de_data %>% filter(comparison == "tissue")

# Top 5 genes per group
top_genes_tissue <- comp_data_tissue %>%
  filter(sig_status != "NS") %>%
  group_by(group) %>%
  slice_min(padj, n = 5, with_ties = FALSE) %>%
  ungroup()

comp_data_tissue <- comp_data_tissue %>%
  mutate(label = ifelse(gene_id %in% top_genes_tissue$gene_id, name, NA_character_))

p_facet_tissue <- ggplot(comp_data_tissue, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(
    data = comp_data_tissue %>% filter(sig_status == "NS"),
    color = "grey80", alpha = 0.4, size = 1
  ) +
  geom_point(
    data = comp_data_tissue %>% filter(sig_status != "NS"),
    aes(color = family),
    alpha = 0.8, size = 2
  ) +
  geom_text_repel(
    data = comp_data_tissue %>% filter(!is.na(label)),
    aes(label = label),
    size = 2.5,
    max.overlaps = 15,
    segment.color = "grey50"
  ) +
  geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = family_colors, name = "Gene Family") +
  facet_wrap(~ group, scales = "free") +
  labs(
    title = "Target Gene Families — Antennae vs Thorax",
    subtitle = "Faceted by functional group | Top 5 genes labeled per group",
    x = expression(log[2]~"Fold Change"),
    y = expression(-log[10]~"adjusted p-value")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 2))

ggsave(file.path(PLOTS_DIR, "volcano_faceted_antennae_vs_thorax.pdf"), 
       p_facet_tissue, width = 14, height = 6)
message("  Saved: volcano_faceted_antennae_vs_thorax.pdf")

# 10) Extract top DE genes for identification.
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 7: Extracting Top DE Genes for Identification")
message(paste(rep("=", 70), collapse = ""))

message("\n[Step 8] Creating table of top DE genes per family...")

# Get top 5 DE genes per family for key comparisons
top_de_for_identification <- all_de_data %>%
  filter(comparison %in% key_comparisons & sig_status != "NS") %>%
  group_by(comparison, group, family) %>%
  slice_min(padj, n = 5, with_ties = FALSE) %>%
  ungroup() %>%
  select(
    comparison, group, family, gene_id, name, product,
    log2FoldChange, padj, sig_status
  ) %>%
  arrange(comparison, group, family, padj)

write_tsv(top_de_for_identification, file.path(OUTPUT_DIR, "top_DE_genes_target_families.tsv"))
message("  Saved: top_DE_genes_target_families.tsv")

# Identify hypothetical proteins among top DE genes
hypothetical_top_de <- top_de_for_identification %>%
  filter(product == "hypothetical protein" | is.na(product) | product == "")

write_tsv(hypothetical_top_de, file.path(OUTPUT_DIR, "hypothetical_proteins_top_DE.tsv"))
message("  Saved: hypothetical_proteins_top_DE.tsv")
message("  Hypothetical proteins among top DE: ", nrow(hypothetical_top_de))

# 11) Summary.
message("\n", paste(rep("=", 70), collapse = ""))
message("SUMMARY")
message(paste(rep("=", 70), collapse = ""))

# Count DE genes per group
summary_stats <- all_de_data %>%
  filter(comparison == "larva_vs_adults" & sig_status != "NS") %>%
  group_by(group) %>%
  summarise(
    n_DE = n(),
    n_up = sum(sig_status == "Up"),
    n_down = sum(sig_status == "Down"),
    pct_up = round(100 * n_up / n_DE, 1),
    .groups = "drop"
  )

message("\n  Larva vs Adults - DE genes by functional group:")
print(summary_stats)

# Hypothetical proteins summary
hypo_summary <- all_de_data %>%
  filter(comparison == "larva_vs_adults" & sig_status != "NS") %>%
  mutate(is_hypothetical = product == "hypothetical protein" | is.na(product)) %>%
  group_by(group) %>%
  summarise(
    n_DE = n(),
    n_hypothetical = sum(is_hypothetical),
    pct_hypothetical = round(100 * n_hypothetical / n_DE, 1),
    .groups = "drop"
  )

message("\n  Hypothetical proteins among DE genes:")
print(hypo_summary)

# 12) Output files.
message("\n", paste(rep("=", 70), collapse = ""))
message("OUTPUT FILES")
message(paste(rep("=", 70), collapse = ""))

message("\nPlots in: ", PLOTS_DIR)
message("  - volcano_Chemosensory_*.pdf (per comparison)")
message("  - volcano_Detoxification_*.pdf (per comparison)")
message("  - volcano_Digestion_*.pdf (per comparison)")
message("  - volcano_all_families_*.pdf (combined)")
message("  - volcano_faceted_larva_vs_adults.pdf")
message("  - volcano_faceted_antennae_vs_thorax.pdf")

message("\nData tables in: ", OUTPUT_DIR)
message("  - top_DE_genes_target_families.tsv (all top DE genes)")
message("  - hypothetical_proteins_top_DE.tsv (for identification)")

message("\n", paste(rep("=", 70), collapse = ""))
message("Next step: Use hypothetical_proteins_top_DE.tsv for BLAST identification")
message(paste(rep("=", 70), collapse = ""))
