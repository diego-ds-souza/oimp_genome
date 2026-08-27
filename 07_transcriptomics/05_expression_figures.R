#!/usr/bin/env Rscript
#
# Expression figures for the host-interface gene families: per-family barplots,
# boxplots and dot plots, volcano and MA plots per contrast, a family-by-
# contrast summary heatmap, lollipop and stacked-bar counts of DE genes, an
# UpSet plot of contrast overlap, PCA of the samples, and per-family heatmaps
# of individual genes.
#
# Family-level expression is the mean normalised expression across the genes of
# a family, which represents typical per-gene deployment and keeps families of
# different sizes comparable.
#
# Run from the repository root:
#   Rscript 07_transcriptomics/05_expression_figures.R

# ------------------------------- settings ------------------------------------
DESEQ_RESULTS <- Sys.getenv("DESEQ_RESULTS",
                            "results/07_transcriptomics/deseq2/deseq2_results.RData")
FAMILIES_DIR  <- Sys.getenv("FAMILIES_DIR", "results/07_transcriptomics/families")
FAMILIES_FILE <- file.path(FAMILIES_DIR, "gene_family_assignments.tsv")
PLOTS_DIR     <- Sys.getenv("PLOTS_DIR", "results/07_transcriptomics/plots")
# -----------------------------------------------------------------------------

dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

# 1) Load packages.
message("\n", paste(rep("=", 70), collapse = ""))
message("Alternative Visualizations Suite")
message(paste(rep("=", 70), collapse = ""))
message("\nLoading packages...")

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
})

# Check for optional packages
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
has_UpSetR <- requireNamespace("UpSetR", quietly = TRUE)
has_ggforce <- requireNamespace("ggforce", quietly = TRUE)

if (has_ggrepel) library(ggrepel)
if (has_patchwork) library(patchwork)
if (has_UpSetR) library(UpSetR)
if (has_ggforce) library(ggforce)

message("  Optional packages: ggrepel=", has_ggrepel, 
        ", patchwork=", has_patchwork, ", UpSetR=", has_UpSetR,
        ", ggforce=", has_ggforce)

# 2) Load data.
message("\n[Step 1] Loading data...")

load(DESEQ_RESULTS)
families_df <- read_tsv(FAMILIES_FILE, show_col_types = FALSE)

message("  Genes: ", nrow(norm_counts_vst))
message("  Samples: ", ncol(norm_counts_vst))
message("  Comparisons: ", length(all_results))

# 3) Define gene families.
FAMILY_GROUPS <- list(
  Chemosensory = c("OBP", "CSP", "OR", "GR", "IR", "SNMP"),
  Detoxification = c("P450", "GST", "UGT", "CCE", "ABC"),
  Digestion = c("SerPro", "CysPro", "AspPro", "MetPro", "CAZyme"),
  Lipid_metabolism = c("FAS", "ELO", "Desat", "ACSL"),
  Terpenoid = c("FPPS", "IspS")
)

# Function to get genes for a family
get_family_genes <- function(family_name, families_dir) {
  fam_file <- file.path(families_dir, paste0("genes_", family_name, ".txt"))
  if (file.exists(fam_file)) {
    genes <- readLines(fam_file)
    return(genes[genes != ""])
  }
  return(character(0))
}

# Collect all genes per family
family_genes <- list()
for (group_name in names(FAMILY_GROUPS)) {
  for (fam in FAMILY_GROUPS[[group_name]]) {
    genes <- get_family_genes(fam, FAMILIES_DIR)
    if (length(genes) > 0) {
      family_genes[[fam]] <- list(genes = genes, group = group_name)
    }
  }
}

message("  Families loaded: ", length(family_genes))

# 4) Prepare sample metadata.
sample_meta <- data.frame(
  sample_id = rownames(sample_info),
  stage = sample_info$stage,
  sex = sample_info$sex,
  tissue = sample_info$structure
) %>%
  mutate(
    condition = paste(stage, sex, sep = "_"),
    full_condition = paste(stage, sex, tissue, sep = "_")
  )

# Color palettes
colors_stage <- c(larva = "#E69F00", adult = "#56B4E9")
colors_sex <- c(male = "#0072B2", female = "#D55E00", unknown = "#999999")
colors_tissue <- c(antenna = "#009E73", thorax = "#CC79A7", head = "#F0E442", midgut = "#E69F00")

# 5) Gene family expression plots.
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 1: Gene Family Expression Plots")
message(paste(rep("=", 70), collapse = ""))

# 6) Calculate family-level expression statistics.
message("\n[Step 2] Calculating family expression statistics...")

calc_family_stats <- function(genes, expr_matrix, sample_meta) {
  genes_present <- genes[genes %in% rownames(expr_matrix)]
  if (length(genes_present) == 0) return(NULL)
  
  mat <- expr_matrix[genes_present, , drop = FALSE]
  
  # Mean expression per sample (across all genes in family)
  sample_means <- colMeans(mat)
  
  result <- data.frame(
    sample_id = names(sample_means),
    mean_expr = sample_means
  ) %>%
    left_join(sample_meta, by = "sample_id")
  
  return(result)
}

# Calculate for all families
family_expr_list <- list()

for (fam_name in names(family_genes)) {
  genes <- family_genes[[fam_name]]$genes
  group <- family_genes[[fam_name]]$group
  
  expr_data <- calc_family_stats(genes, norm_counts_vst, sample_meta)
  
  if (!is.null(expr_data)) {
    expr_data$family <- fam_name
    expr_data$group <- group
    expr_data$n_genes <- sum(genes %in% rownames(norm_counts_vst))
    family_expr_list[[fam_name]] <- expr_data
  }
}

family_expr_df <- bind_rows(family_expr_list)
message("  Families with expression data: ", n_distinct(family_expr_df$family))

# Calculate summary statistics per family/condition
family_summary <- family_expr_df %>%
  group_by(family, group, stage, sex, tissue) %>%
  summarise(
    mean = mean(mean_expr),
    sd = sd(mean_expr),
    se = sd(mean_expr) / sqrt(n()),
    n = n(),
    n_genes = first(n_genes),
    .groups = "drop"
  ) %>%
  mutate(
    se = ifelse(is.na(se), 0, se),  # Handle n=1 cases
    condition = paste(stage, sex, sep = "_")
  )

# Also by stage only
family_summary_stage <- family_expr_df %>%
  group_by(family, group, stage) %>%
  summarise(
    mean = mean(mean_expr),
    sd = sd(mean_expr),
    se = sd(mean_expr) / sqrt(n()),
    n = n(),
    n_genes = first(n_genes),
    .groups = "drop"
  ) %>%
  mutate(se = ifelse(is.na(se), 0, se))

# By sex (adults only)
family_summary_sex <- family_expr_df %>%
  filter(stage == "adult") %>%
  group_by(family, group, sex) %>%
  summarise(
    mean = mean(mean_expr),
    sd = sd(mean_expr),
    se = sd(mean_expr) / sqrt(n()),
    n = n(),
    n_genes = first(n_genes),
    .groups = "drop"
  ) %>%
  mutate(se = ifelse(is.na(se), 0, se))

# By tissue (adults only)
family_summary_tissue <- family_expr_df %>%
  filter(stage == "adult") %>%
  group_by(family, group, tissue) %>%
  summarise(
    mean = mean(mean_expr),
    sd = sd(mean_expr),
    se = sd(mean_expr) / sqrt(n()),
    n = n(),
    n_genes = first(n_genes),
    .groups = "drop"
  ) %>%
  mutate(se = ifelse(is.na(se), 0, se))

# 7) Barplots - mean expression per family.
message("\n[Step 3] Creating barplots...")

# Barplot by stage
p_bar_stage <- ggplot(family_summary_stage, 
                       aes(x = reorder(family, -mean), y = mean, fill = stage)) +

geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(0.8), width = 0.25) +
  facet_wrap(~group, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = colors_stage) +
  labs(
    title = "Gene Family Expression by Developmental Stage (Barplot)",
    subtitle = "Mean VST-normalized expression ± SE",
    x = "Gene Family",
    y = "Mean Expression",
    fill = "Stage"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(PLOTS_DIR, "01_barplot_expression_by_stage.pdf"), 
       p_bar_stage, width = 14, height = 6)
message("  Saved: 01_barplot_expression_by_stage.pdf")

# Barplot by sex (adults)
p_bar_sex <- ggplot(family_summary_sex, 
                     aes(x = reorder(family, -mean), y = mean, fill = sex)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(0.8), width = 0.25) +
  facet_wrap(~group, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = colors_sex) +
  labs(
    title = "Gene Family Expression by Sex (Barplot, Adults Only)",
    subtitle = "Mean VST-normalized expression ± SE",
    x = "Gene Family",
    y = "Mean Expression",
    fill = "Sex"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(PLOTS_DIR, "02_barplot_expression_by_sex.pdf"), 
       p_bar_sex, width = 14, height = 6)
message("  Saved: 02_barplot_expression_by_sex.pdf")

# Barplot by tissue (adults)
p_bar_tissue <- ggplot(family_summary_tissue, 
                        aes(x = reorder(family, -mean), y = mean, fill = tissue)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(0.8), width = 0.25) +
  facet_wrap(~group, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = colors_tissue) +
  labs(
    title = "Gene Family Expression by Tissue (Barplot, Adults Only)",
    subtitle = "Mean VST-normalized expression ± SE",
    x = "Gene Family",
    y = "Mean Expression",
    fill = "Tissue"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(PLOTS_DIR, "03_barplot_expression_by_tissue.pdf"), 
       p_bar_tissue, width = 14, height = 6)
message("  Saved: 03_barplot_expression_by_tissue.pdf")

# 8) Boxplots - distribution of expression.
message("\n[Step 4] Creating boxplots...")

# Boxplot by stage
p_box_stage <- ggplot(family_expr_df, 
                       aes(x = reorder(family, -mean_expr), y = mean_expr, fill = stage)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.size = 1) +
  facet_wrap(~group, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = colors_stage) +
  labs(
    title = "Gene Family Expression by Developmental Stage (Boxplot)",
    subtitle = "Distribution of mean family expression across samples",
    x = "Gene Family",
    y = "Mean Expression (VST)",
    fill = "Stage"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(PLOTS_DIR, "04_boxplot_expression_by_stage.pdf"), 
       p_box_stage, width = 14, height = 6)
message("  Saved: 04_boxplot_expression_by_stage.pdf")

# Boxplot by sex (adults)
p_box_sex <- ggplot(family_expr_df %>% filter(stage == "adult"), 
                     aes(x = reorder(family, -mean_expr), y = mean_expr, fill = sex)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.size = 1) +
  facet_wrap(~group, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = colors_sex) +
  labs(
    title = "Gene Family Expression by Sex (Boxplot, Adults Only)",
    subtitle = "Distribution of mean family expression across samples",
    x = "Gene Family",
    y = "Mean Expression (VST)",
    fill = "Sex"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(PLOTS_DIR, "05_boxplot_expression_by_sex.pdf"), 
       p_box_sex, width = 14, height = 6)
message("  Saved: 05_boxplot_expression_by_sex.pdf")

# Boxplot by tissue (adults)
p_box_tissue <- ggplot(family_expr_df %>% filter(stage == "adult"), 
                        aes(x = reorder(family, -mean_expr), y = mean_expr, fill = tissue)) +
  geom_boxplot(position = position_dodge(0.8), width = 0.7, outlier.size = 1) +
  facet_wrap(~group, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = colors_tissue) +
  labs(
    title = "Gene Family Expression by Tissue (Boxplot, Adults Only)",
    subtitle = "Distribution of mean family expression across samples",
    x = "Gene Family",
    y = "Mean Expression (VST)",
    fill = "Tissue"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(PLOTS_DIR, "06_boxplot_expression_by_tissue.pdf"), 
       p_box_tissue, width = 14, height = 6)
message("  Saved: 06_boxplot_expression_by_tissue.pdf")

# 9) Dot plots with error bars - mean ± SE (better for small n).
message("\n[Step 5] Creating dot plots with error bars...")

# Dot plot by stage
p_dot_stage <- ggplot(family_summary_stage, 
                       aes(x = reorder(family, -mean), y = mean, color = stage)) +
  geom_point(position = position_dodge(0.6), size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(0.6), width = 0.3) +
  facet_wrap(~group, scales = "free_x", nrow = 1) +
  scale_color_manual(values = colors_stage) +
  labs(
    title = "Gene Family Expression by Developmental Stage (Dot Plot)",
    subtitle = "Mean ± SE; better visualization for small sample sizes",
    x = "Gene Family",
    y = "Mean Expression (VST)",
    color = "Stage"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(PLOTS_DIR, "07_dotplot_expression_by_stage.pdf"), 
       p_dot_stage, width = 14, height = 6)
message("  Saved: 07_dotplot_expression_by_stage.pdf")

# Dot plot by sex
p_dot_sex <- ggplot(family_summary_sex, 
                     aes(x = reorder(family, -mean), y = mean, color = sex)) +
  geom_point(position = position_dodge(0.6), size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(0.6), width = 0.3) +
  facet_wrap(~group, scales = "free_x", nrow = 1) +
  scale_color_manual(values = colors_sex) +
  labs(
    title = "Gene Family Expression by Sex (Dot Plot, Adults Only)",
    subtitle = "Mean ± SE",
    x = "Gene Family",
    y = "Mean Expression (VST)",
    color = "Sex"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(PLOTS_DIR, "08_dotplot_expression_by_sex.pdf"), 
       p_dot_sex, width = 14, height = 6)
message("  Saved: 08_dotplot_expression_by_sex.pdf")

# Dot plot by tissue
p_dot_tissue <- ggplot(family_summary_tissue, 
                        aes(x = reorder(family, -mean), y = mean, color = tissue)) +
  geom_point(position = position_dodge(0.6), size = 3) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(0.6), width = 0.3) +
  facet_wrap(~group, scales = "free_x", nrow = 1) +
  scale_color_manual(values = colors_tissue) +
  labs(
    title = "Gene Family Expression by Tissue (Dot Plot, Adults Only)",
    subtitle = "Mean ± SE",
    x = "Gene Family",
    y = "Mean Expression (VST)",
    color = "Tissue"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(PLOTS_DIR, "09_dotplot_expression_by_tissue.pdf"), 
       p_dot_tissue, width = 14, height = 6)
message("  Saved: 09_dotplot_expression_by_tissue.pdf")

# 10) Differential expression plots.
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 2: Differential Expression Plots")
message(paste(rep("=", 70), collapse = ""))

# 11) Volcano plots.
message("\n[Step 6] Creating volcano plots...")

create_volcano_plot <- function(de_result, comparison_name, 
                                 padj_thresh = 0.05, lfc_thresh = 1,
                                 families_df = NULL) {
  
  de_df <- as.data.frame(de_result)
  de_df$gene_id <- rownames(de_df)
  
  # Add significance category
  de_df <- de_df %>%
    mutate(
      significance = case_when(
        is.na(padj) ~ "NS",
        padj >= padj_thresh ~ "NS",
        abs(log2FoldChange) < lfc_thresh ~ "Sig (small FC)",
        log2FoldChange >= lfc_thresh ~ "Up",
        log2FoldChange <= -lfc_thresh ~ "Down",
        TRUE ~ "NS"
      ),
      neg_log10_padj = -log10(pvalue)
    )
  
  # Add family info if available
  if (!is.null(families_df)) {
    de_df <- de_df %>%
      left_join(families_df %>% select(gene_id, families), by = "gene_id")
  }
  
  # Cap extreme values for visualization
  de_df$neg_log10_padj[de_df$neg_log10_padj > 50] <- 50
  
  # Colors
  sig_colors <- c("Up" = "#D55E00", "Down" = "#0072B2", 
                  "Sig (small FC)" = "#999999", "NS" = "#CCCCCC")
  
  # Count significant genes
  n_up <- sum(de_df$significance == "Up", na.rm = TRUE)
  n_down <- sum(de_df$significance == "Down", na.rm = TRUE)
  
  p <- ggplot(de_df, aes(x = log2FoldChange, y = neg_log10_padj, color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", color = "grey40") +
    scale_color_manual(values = sig_colors) +
    labs(
      title = paste0("Volcano Plot: ", comparison_name),
      subtitle = paste0("Up: ", n_up, " | Down: ", n_down, 
                        " (padj < ", padj_thresh, ", |log2FC| > ", lfc_thresh, ")"),
      x = "log2 Fold Change",
      y = "-log10(p-value)",
      color = "Significance"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    ) +
    xlim(c(-max(abs(de_df$log2FoldChange), na.rm = TRUE) - 0.5,
           max(abs(de_df$log2FoldChange), na.rm = TRUE) + 0.5))
  
  # Add labels for top genes if ggrepel available
  if (has_ggrepel) {
    top_genes <- de_df %>%
      filter(significance %in% c("Up", "Down")) %>%
      arrange(padj) %>%
      head(10)
    
    if (nrow(top_genes) > 0) {
      p <- p + geom_text_repel(
        data = top_genes,
        aes(label = gene_id),
        size = 2.5,
        max.overlaps = 15,
        segment.color = "grey50"
      )
    }
  }
  
  return(p)
}

# Create volcano plots for each comparison
comparison_labels <- c(
  "antennae_MvF" = "Male vs Female (Antennae)",
  "thorax_MvF" = "Male vs Female (Thorax)",
  "adults_MvF" = "Male vs Female (All Adults)",
  "tissue" = "Antennae vs Thorax",
  "larva_vs_male" = "Larva vs Adult Male [DESCRIPTIVE]",
  "larva_vs_female" = "Larva vs Adult Female [DESCRIPTIVE]",
  "larva_vs_adults" = "Larva vs All Adults [DESCRIPTIVE]"
)

for (comp_name in names(all_results)) {
  label <- ifelse(comp_name %in% names(comparison_labels), 
                  comparison_labels[comp_name], comp_name)
  
  p <- create_volcano_plot(
    all_results[[comp_name]]$results,
    label,
    families_df = families_df
  )
  
  filename <- paste0("10_volcano_", comp_name, ".pdf")
  ggsave(file.path(PLOTS_DIR, filename), p, width = 10, height = 8)
  message("  Saved: ", filename)
}

# -----------------------------------------------------------------------------
# 2B. DOT PLOT FOR DE RESULTS (Family × Comparison)
# UPDATED: Stronger orange color that doesn't fade to white
# -----------------------------------------------------------------------------

message("\n[Step 7] Creating DE summary dot plot...")

# Calculate DE stats per family per comparison
de_family_stats <- list()

for (comp_name in names(all_results)) {
  de_res <- all_results[[comp_name]]$results
  
  for (fam_name in names(family_genes)) {
    genes <- family_genes[[fam_name]]$genes
    genes_present <- genes[genes %in% rownames(de_res)]
    
    if (length(genes_present) == 0) next
    
    de_sub <- de_res[genes_present, ]
    
    sig <- !is.na(de_sub$padj) & de_sub$padj < 0.05
    up <- sig & de_sub$log2FoldChange > 1
    down <- sig & de_sub$log2FoldChange < -1
    
    mean_lfc <- mean(de_sub$log2FoldChange, na.rm = TRUE)
    
    de_family_stats[[length(de_family_stats) + 1]] <- data.frame(
      comparison = comp_name,
      family = fam_name,
      group = family_genes[[fam_name]]$group,
      total = length(genes_present),
      significant = sum(sig),
      up = sum(up),
      down = sum(down),
      pct_sig = 100 * sum(sig) / length(genes_present),
      mean_log2FC = mean_lfc
    )
  }
}

de_family_df <- bind_rows(de_family_stats)

# Create readable comparison labels
de_family_df <- de_family_df %>%
  mutate(
    comparison_label = case_when(
      comparison == "antennae_MvF" ~ "M vs F\n(antennae)",
      comparison == "thorax_MvF" ~ "M vs F\n(thorax)",
      comparison == "adults_MvF" ~ "M vs F\n(pooled)",
      comparison == "tissue" ~ "Ant vs Th",
      comparison == "larva_vs_male" ~ "L vs M*",
      comparison == "larva_vs_female" ~ "L vs F*",
      comparison == "larva_vs_adults" ~ "L vs A*",
      TRUE ~ comparison
    )
  )

# Dot plot: size = # significant, color = mean log2FC
# UPDATED COLOR SCALE: Using a diverging palette with stronger colors
# Changed from gradient2 with white midpoint to a custom palette with light gray midpoint
p_de_dot <- ggplot(de_family_df %>% filter(significant > 0), 
                    aes(x = comparison_label, y = family)) +
  geom_point(aes(size = significant, color = mean_log2FC)) +
  scale_size_continuous(range = c(2, 10), name = "# Significant\nGenes") +
  # UPDATED: Stronger color scale - orange (#E65100) to light gray (#E0E0E0) to blue (#1565C0)
  scale_color_gradient2(
    low = "#1565C0",      # Strong blue for negative (down in first condition)
    mid = "#E0E0E0",      # Light gray instead of white for better visibility
    high = "#DE2D26",     # Strong orange for positive (up in first condition)
    midpoint = 0, 
    name = "Mean\nlog2FC",
    limits = c(min(de_family_df$mean_log2FC, na.rm = TRUE),
               max(de_family_df$mean_log2FC, na.rm = TRUE))
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Differential Expression by Gene Family",
    subtitle = "Size = number of significant genes; Color = direction of change\n* = descriptive only (limited replicates)",
    x = "Comparison",
    y = "Gene Family"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
    strip.background = element_rect(fill = "grey90"),
    strip.text.y = element_text(angle = 0),
    plot.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey90")
  )

ggsave(file.path(PLOTS_DIR, "11_dotplot_DE_family_comparison.pdf"), 
       p_de_dot, width = 10, height = 10)
message("  Saved: 11_dotplot_DE_family_comparison.pdf")

# 12) Table visualizations.
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 3: Table Visualizations")
message(paste(rep("=", 70), collapse = ""))

# 13) Summary heatmap (families × comparisons).
message("\n[Step 8] Creating summary heatmap...")

# Pivot to matrix format
heatmap_matrix <- de_family_df %>%
  select(family, comparison, significant) %>%
  pivot_wider(names_from = comparison, values_from = significant, values_fill = 0) %>%
  column_to_rownames("family") %>%
  as.matrix()

# Reorder columns
col_order <- c("antennae_MvF", "thorax_MvF", "adults_MvF", "tissue", 
               "larva_vs_male", "larva_vs_female", "larva_vs_adults")
col_order <- col_order[col_order %in% colnames(heatmap_matrix)]
heatmap_matrix <- heatmap_matrix[, col_order]

# Create annotation for rows (gene family groups)
row_annot <- de_family_df %>%
  select(family, group) %>%
  distinct() %>%
  column_to_rownames("family")

# Colors for annotation
group_colors <- c(
  Chemosensory = "#009E73",
  Detoxification = "#D55E00",
  Digestion = "#0072B2",
  Lipid_metabolism = "#CC79A7",
  Terpenoid = "#E69F00"
)

annot_colors_hm <- list(group = group_colors)

# Column labels
col_labels <- c(
  "antennae_MvF" = "M vs F (ant)",
  "thorax_MvF" = "M vs F (th)",
  "adults_MvF" = "M vs F (all)",
  "tissue" = "Ant vs Th",
  "larva_vs_male" = "L vs M*",
  "larva_vs_female" = "L vs F*",
  "larva_vs_adults" = "L vs A*"
)

colnames(heatmap_matrix) <- col_labels[colnames(heatmap_matrix)]

# Create heatmap
pdf(file.path(PLOTS_DIR, "12_heatmap_DE_summary.pdf"), width = 10, height = 10)
pheatmap(
  heatmap_matrix,
  main = "Number of Significant DE Genes per Family\n(padj < 0.05, |log2FC| > 1)",
  color = colorRampPalette(c("white", "#FEE0D2", "#FC9272", "#DE2D26"))(50),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_row = row_annot,
  annotation_colors = annot_colors_hm,
  display_numbers = TRUE,
  number_format = "%d",
  fontsize = 10,
  fontsize_number = 8,
  border_color = "grey80"
)
dev.off()
message("  Saved: 12_heatmap_DE_summary.pdf")

# 14) Lollipop chart (DE gene counts by family).
message("\n[Step 9] Creating lollipop charts...")

# For adults_MvF comparison
lollipop_data <- de_family_df %>%
  filter(comparison == "adults_MvF") %>%
  arrange(desc(significant))

p_lollipop <- ggplot(lollipop_data, aes(x = reorder(family, significant), y = significant)) +
  geom_segment(aes(xend = family, yend = 0, color = group), size = 1.5) +
  geom_point(aes(color = group), size = 4) +
  geom_text(aes(label = significant), hjust = -0.5, size = 3) +
  coord_flip() +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Significant DE Genes by Family",
    subtitle = "Comparison: Male vs Female (All Adults)",
    x = "Gene Family",
    y = "Number of Significant Genes",
    color = "Functional Group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  ylim(0, max(lollipop_data$significant) * 1.1)

ggsave(file.path(PLOTS_DIR, "13_lollipop_DE_adults_MvF.pdf"), 
       p_lollipop, width = 8, height = 10)
message("  Saved: 13_lollipop_DE_adults_MvF.pdf")

# For larva_vs_adults comparison
lollipop_data_larva <- de_family_df %>%
  filter(comparison == "larva_vs_adults") %>%
  arrange(desc(significant))

p_lollipop_larva <- ggplot(lollipop_data_larva, 
                            aes(x = reorder(family, significant), y = significant)) +
  geom_segment(aes(xend = family, yend = 0, color = group), size = 1.5) +
  geom_point(aes(color = group), size = 4) +
  geom_text(aes(label = significant), hjust = -0.5, size = 3) +
  coord_flip() +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Significant DE Genes by Family [DESCRIPTIVE]",
    subtitle = "Comparison: Larva vs All Adults (interpret with caution)",
    x = "Gene Family",
    y = "Number of Significant Genes",
    color = "Functional Group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  ylim(0, max(lollipop_data_larva$significant) * 1.1)

ggsave(file.path(PLOTS_DIR, "14_lollipop_DE_larva_vs_adults.pdf"), 
       p_lollipop_larva, width = 8, height = 10)
message("  Saved: 14_lollipop_DE_larva_vs_adults.pdf")

# 15) Stacked bar chart (up vs down per family).
message("\n[Step 10] Creating stacked bar charts...")

# Reshape data for stacked bar
stacked_data <- de_family_df %>%
  select(comparison, family, group, up, down) %>%
  pivot_longer(cols = c(up, down), names_to = "direction", values_to = "count") %>%
  mutate(
    direction = factor(direction, levels = c("up", "down")),
    count_signed = ifelse(direction == "down", -count, count)
  )

# For adults_MvF
p_stacked <- ggplot(stacked_data %>% filter(comparison == "adults_MvF"),
                     aes(x = reorder(family, abs(count_signed)), y = count_signed, fill = direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("up" = "#D55E00", "down" = "#0072B2"),
                    labels = c("up" = "Upregulated", "down" = "Downregulated")) +
  geom_hline(yintercept = 0, color = "black") +
  labs(
    title = "Direction of DE by Gene Family",
    subtitle = "Comparison: Male vs Female (All Adults)",
    x = "Gene Family",
    y = "Number of Genes (negative = downregulated)",
    fill = "Direction"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(PLOTS_DIR, "15_stacked_bar_DE_direction_adults_MvF.pdf"), 
       p_stacked, width = 8, height = 10)
message("  Saved: 15_stacked_bar_DE_direction_adults_MvF.pdf")

# For larva_vs_adults
p_stacked_larva <- ggplot(stacked_data %>% filter(comparison == "larva_vs_adults"),
                           aes(x = reorder(family, abs(count_signed)), y = count_signed, fill = direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("up" = "#D55E00", "down" = "#0072B2"),
                    labels = c("up" = "Up in larvae", "down" = "Down in larvae")) +
  geom_hline(yintercept = 0, color = "black") +
  labs(
    title = "Direction of DE by Gene Family [DESCRIPTIVE]",
    subtitle = "Comparison: Larva vs All Adults (up = higher in larvae)",
    x = "Gene Family",
    y = "Number of Genes",
    fill = "Direction"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(PLOTS_DIR, "16_stacked_bar_DE_direction_larva_vs_adults.pdf"), 
       p_stacked_larva, width = 8, height = 10)
message("  Saved: 16_stacked_bar_DE_direction_larva_vs_adults.pdf")

# 16) UpSet plot (overlap between comparisons).
message("\n[Step 11] Creating UpSet plot...")

if (has_UpSetR) {
  # Get significant genes for each comparison
  sig_genes_list <- list()
  
  for (comp_name in names(all_results)) {
    de_res <- all_results[[comp_name]]$results
    sig_genes <- rownames(de_res)[!is.na(de_res$padj) & de_res$padj < 0.05]
    
    # Use readable names
    label <- comparison_labels[comp_name]
    if (is.na(label)) label <- comp_name
    
    sig_genes_list[[label]] <- sig_genes
  }
  
  # Create UpSet plot
  pdf(file.path(PLOTS_DIR, "17_upset_DE_overlap.pdf"), width = 12, height = 8)
  print(upset(
    fromList(sig_genes_list),
    nsets = length(sig_genes_list),
    order.by = "freq",
    decreasing = TRUE,
    mainbar.y.label = "Number of Shared DE Genes",
    sets.x.label = "DE Genes per Comparison",
    text.scale = 1.2,
    point.size = 3,
    line.size = 1
  ))
  dev.off()
  message("  Saved: 17_upset_DE_overlap.pdf")
} else {
  message("  Skipping UpSet plot (UpSetR package not installed)")
}

# 17) Sample relationship plots.
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 4: Sample Relationship Plots")
message(paste(rep("=", 70), collapse = ""))

# 18) PCA plots.
message("\n[Step 12] Creating PCA plots...")

pca_data <- prcomp(t(norm_counts_vst), scale. = FALSE)
var_explained <- round(100 * pca_data$sdev^2 / sum(pca_data$sdev^2), 1)

pca_df <- data.frame(
  PC1 = pca_data$x[, 1],
  PC2 = pca_data$x[, 2],
  PC3 = pca_data$x[, 3],
  sample_id = rownames(pca_data$x)
) %>%
  left_join(sample_meta, by = "sample_id")

# All samples, sex legend without NA.
# Filter out samples where sex is NA or "unknown" for the shape aesthetic
# But keep all samples for plotting
pca_df_clean <- pca_df %>%
  mutate(
    # Replace NA and "unknown" with a specific value that we'll handle
    sex_clean = case_when(
      is.na(sex) ~ "not determined",
      sex == "unknown" ~ "not determined",
      TRUE ~ sex
    )
  )

# Check if we have any "not determined" - if all samples have sex, we can exclude from legend
# If larvae have unknown sex, we need to decide how to handle
n_undetermined <- sum(pca_df_clean$sex_clean == "not determined")

if (n_undetermined > 0) {
  # Option: Don't map sex to shape for samples without sex info
  # Create plot with only stage colored, and sex only for adults
  
  # For larvae, we won't show sex shape
  pca_df_plot <- pca_df_clean %>%
    mutate(
      # Only show sex shapes for adults
      sex_for_plot = ifelse(stage == "adult", sex_clean, NA_character_)
    )
  
  # Create base plot - all points colored by stage
  p_pca <- ggplot(pca_df_plot, aes(x = PC1, y = PC2)) +
    # First layer: all points colored by stage
    geom_point(aes(color = stage), size = 5) +
    # Second layer: adult points with sex shapes (overplot)
    geom_point(data = pca_df_plot %>% filter(stage == "adult"),
               aes(color = stage, shape = sex_for_plot), size = 5) +
    scale_color_manual(values = colors_stage) +
    scale_shape_manual(values = c(male = 17, female = 16), na.translate = FALSE) +
    labs(
      title = "PCA of Gene Expression",
      subtitle = "All samples, VST-normalized counts",
      x = paste0("PC1 (", var_explained[1], "% variance)"),
      y = paste0("PC2 (", var_explained[2], "% variance)"),
      color = "Stage",
      shape = "Sex"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    guides(
      color = guide_legend(order = 1),
      shape = guide_legend(order = 2)
    )
} else {
  # All samples have sex info - simple plot
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = stage, shape = sex), size = 5) +
    scale_color_manual(values = colors_stage) +
    scale_shape_manual(values = c(male = 17, female = 16)) +
    labs(
      title = "PCA of Gene Expression",
      subtitle = "All samples, VST-normalized counts",
      x = paste0("PC1 (", var_explained[1], "% variance)"),
      y = paste0("PC2 (", var_explained[2], "% variance)"),
      color = "Stage",
      shape = "Sex"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

if (has_ggrepel) {
  p_pca <- p_pca + geom_text_repel(aes(label = sample_id), size = 3, max.overlaps = 20)
} else {
  p_pca <- p_pca + geom_text(aes(label = sample_id), vjust = -1, size = 3)
}

ggsave(file.path(PLOTS_DIR, "18_PCA_all_samples.pdf"), p_pca, width = 10, height = 8)
message("  Saved: 18_PCA_all_samples.pdf")

# The same, with cluster ellipses.
if (has_ggforce) {
  # Create plot with ellipses around stage clusters
  p_pca_ellipse <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    # Add ellipses first (behind points)
    geom_mark_ellipse(aes(fill = stage, label = stage), 
                      alpha = 0.15, 
                      expand = unit(3, "mm"),
                      label.fontsize = 10,
                      label.fill = "white",
                      con.cap = 0) +
    # Add points
    geom_point(aes(color = stage), size = 5) +
    scale_color_manual(values = colors_stage) +
    scale_fill_manual(values = colors_stage) +
    labs(
      title = "PCA of Gene Expression",
      subtitle = "All samples with cluster ellipses by developmental stage",
      x = paste0("PC1 (", var_explained[1], "% variance)"),
      y = paste0("PC2 (", var_explained[2], "% variance)"),
      color = "Stage"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    guides(fill = "none")  # Hide fill legend (redundant with color)
  
  if (has_ggrepel) {
    p_pca_ellipse <- p_pca_ellipse + 
      geom_text_repel(aes(label = sample_id), size = 3, max.overlaps = 20)
  } else {
    p_pca_ellipse <- p_pca_ellipse + 
      geom_text(aes(label = sample_id), vjust = -1, size = 3)
  }
  
  ggsave(file.path(PLOTS_DIR, "18b_PCA_all_samples_with_ellipses.pdf"), 
         p_pca_ellipse, width = 10, height = 8)
  message("  Saved: 18b_PCA_all_samples_with_ellipses.pdf")
} else {
  message("  Skipping 18b (ggforce not installed for ellipses)")
}

# Adults only, PCA recomputed on that subset.
# Get adult sample IDs
adult_samples <- sample_meta %>% filter(stage == "adult") %>% pull(sample_id)

# Run NEW PCA on adults only
pca_adults_data <- prcomp(t(norm_counts_vst[, adult_samples]), scale. = FALSE)
var_explained_adults <- round(100 * pca_adults_data$sdev^2 / sum(pca_adults_data$sdev^2), 1)

pca_adults_df <- data.frame(
  PC1 = pca_adults_data$x[, 1],
  PC2 = pca_adults_data$x[, 2],
  sample_id = rownames(pca_adults_data$x)
) %>%
  left_join(sample_meta, by = "sample_id")

message("  Adults-only PCA: PC1=", var_explained_adults[1], "%, PC2=", var_explained_adults[2], "%")

p_pca_tissue <- ggplot(pca_adults_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = tissue, shape = sex), size = 5) +
  scale_color_manual(values = colors_tissue) +
  scale_shape_manual(values = c(male = 17, female = 16)) +
  labs(
    title = "PCA of Gene Expression (Adults Only)",
    subtitle = "Colored by tissue",
    x = paste0("PC1 (", var_explained_adults[1], "% variance)"),
    y = paste0("PC2 (", var_explained_adults[2], "% variance)"),
    color = "Tissue",
    shape = "Sex"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

if (has_ggrepel) {
  p_pca_tissue <- p_pca_tissue + geom_text_repel(aes(label = sample_id), size = 3)
}

ggsave(file.path(PLOTS_DIR, "19_PCA_adults_by_tissue.pdf"), p_pca_tissue, width = 10, height = 8)
message("  Saved: 19_PCA_adults_by_tissue.pdf")

# The same, with cluster ellipses.
if (has_ggforce) {
  p_pca_tissue_ellipse <- ggplot(pca_adults_df, aes(x = PC1, y = PC2)) +
    # Add ellipses around tissue clusters
    geom_mark_ellipse(aes(fill = tissue, label = tissue), 
                      alpha = 0.15, 
                      expand = unit(3, "mm"),
                      label.fontsize = 10,
                      label.fill = "white",
                      con.cap = 0) +
    # Add points
    geom_point(aes(color = tissue, shape = sex), size = 5) +
    scale_color_manual(values = colors_tissue) +
    scale_fill_manual(values = colors_tissue) +
    scale_shape_manual(values = c(male = 17, female = 16)) +
    labs(
      title = "PCA of Gene Expression (Adults Only)",
      subtitle = "With cluster ellipses by tissue type",
      x = paste0("PC1 (", var_explained_adults[1], "% variance)"),
      y = paste0("PC2 (", var_explained_adults[2], "% variance)"),
      color = "Tissue",
      shape = "Sex"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    guides(fill = "none")
  
  if (has_ggrepel) {
    p_pca_tissue_ellipse <- p_pca_tissue_ellipse + 
      geom_text_repel(aes(label = sample_id), size = 3)
  }
  
  ggsave(file.path(PLOTS_DIR, "19b_PCA_adults_by_tissue_with_ellipses.pdf"), 
         p_pca_tissue_ellipse, width = 10, height = 8)
  message("  Saved: 19b_PCA_adults_by_tissue_with_ellipses.pdf")
} else {
  message("  Skipping 19b (ggforce not installed for ellipses)")
}

# 19) Gene family heatmaps (individual gene expression).
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 5: Gene Family Heatmaps")
message(paste(rep("=", 70), collapse = ""))

message("\n[Step 13] Creating gene family heatmaps...")

# Create sample annotation for heatmaps
sample_annot <- data.frame(
  row.names = colnames(norm_counts_vst),
  Stage = sample_info$stage,
  Sex = sample_info$sex,
  Tissue = sample_info$structure
)

# Replace NA/unknown in Sex for cleaner annotation
sample_annot$Sex[is.na(sample_annot$Sex) | sample_annot$Sex == "unknown"] <- "ND"

# Define colors for heatmap annotations
annot_colors_hm <- list(
  Stage = c(larva = "#E69F00", adult = "#56B4E9"),
  Sex = c(male = "#0072B2", female = "#D55E00", ND = "#CCCCCC"),
  Tissue = c(antenna = "#009E73", thorax = "#CC79A7", head = "#F0E442", midgut = "#E69F00")
)

# Helper function to create and save heatmap
create_family_heatmap <- function(genes, expr_matrix, family_name, 
                                   sample_annotation, annotation_colors,
                                   output_dir, show_rownames = TRUE) {
  
  # Filter to genes present in expression matrix
  genes_present <- genes[genes %in% rownames(expr_matrix)]
  
  if (length(genes_present) == 0) {
    message("  No genes found for ", family_name)
    return(NULL)
  }
  
  if (length(genes_present) < length(genes)) {
    message("  ", family_name, ": ", length(genes_present), "/", length(genes), 
            " genes found")
  }
  
  # Extract expression data
  mat <- expr_matrix[genes_present, , drop = FALSE]
  
  # Scale by row (z-score) for visualization
  mat_scaled <- t(scale(t(mat)))
  
  # Handle any NA/Inf from scaling
  mat_scaled[is.na(mat_scaled)] <- 0
  mat_scaled[is.infinite(mat_scaled)] <- 0
  
  # Determine row labels
  if (nrow(mat_scaled) > 50) {
    show_rownames <- FALSE
    message("    More than 50 genes, hiding row names")
  }
  
  # Create title
  title <- paste0(family_name, " (n=", nrow(mat_scaled), " genes)")
  
  # Color palette (blue-white-red diverging)
  color_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
  
  # Calculate height based on number of genes
  plot_height <- max(6, nrow(mat_scaled) * 0.15 + 3)
  
  # Create and save heatmap
  filename <- file.path(output_dir, paste0("20_heatmap_", 
                                            gsub(" ", "_", gsub("[^[:alnum:] ]", "", family_name)), 
                                            ".pdf"))
  
  pdf(filename, width = 10, height = plot_height)
  
  pheatmap(
    mat_scaled,
    main = title,
    color = color_palette,
    breaks = seq(-3, 3, length.out = 101),
    cluster_rows = TRUE,
    cluster_cols = FALSE,  # Keep sample order
    show_rownames = show_rownames,
    show_colnames = TRUE,
    annotation_col = sample_annotation,
    annotation_colors = annotation_colors,
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 10,
    border_color = NA
  )
  
  dev.off()
  message("  Saved: ", basename(filename))
  
  return(list(genes = genes_present, n_genes = length(genes_present)))
}

# Priority families to create heatmaps for (most biologically relevant)
priority_families <- c("OBP", "OR", "CSP", "P450", "GST", "UGT", "SerPro", "CAZyme")

# Create heatmaps for priority families
for (fam_name in priority_families) {
  if (fam_name %in% names(family_genes)) {
    genes <- family_genes[[fam_name]]$genes
    group <- family_genes[[fam_name]]$group
    
    create_family_heatmap(
      genes = genes,
      expr_matrix = norm_counts_vst,
      family_name = paste0(group, " - ", fam_name),
      sample_annotation = sample_annot,
      annotation_colors = annot_colors_hm,
      output_dir = PLOTS_DIR
    )
  }
}

# Also create combined heatmaps for each major functional group
message("\n  Creating combined group heatmaps...")

for (group_name in c("Chemosensory", "Detoxification", "Digestion")) {
  # Collect all genes from families in this group
  group_families <- FAMILY_GROUPS[[group_name]]
  all_group_genes <- c()
  
  for (fam in group_families) {
    if (fam %in% names(family_genes)) {
      all_group_genes <- c(all_group_genes, family_genes[[fam]]$genes)
    }
  }
  
  if (length(all_group_genes) > 0) {
    create_family_heatmap(
      genes = all_group_genes,
      expr_matrix = norm_counts_vst,
      family_name = paste0(group_name, " - All Families"),
      sample_annotation = sample_annot,
      annotation_colors = annot_colors_hm,
      output_dir = PLOTS_DIR
    )
  }
}

# 20) Save data and summary.
message("\n[Step 14] Saving data objects...")

save(
  family_expr_df,
  family_summary,
  family_summary_stage,
  family_summary_sex,
  family_summary_tissue,
  de_family_df,
  pca_df,
  file = file.path(PLOTS_DIR, "visualization_data.RData")
)
message("  Saved: visualization_data.RData")

# Write summary tables
write_csv(de_family_df, file.path(PLOTS_DIR, "DE_family_summary_table.csv"))
message("  Saved: DE_family_summary_table.csv")

# 21) Final summary.
message("\n", paste(rep("=", 70), collapse = ""))
message("VISUALIZATION COMPLETE")
message(paste(rep("=", 70), collapse = ""))

message("\nAll plots saved to: ", PLOTS_DIR)
message("\nPlots generated:")
message("  EXPRESSION PLOTS:")
message("    01-03: Barplots (stage, sex, tissue)")
message("    04-06: Boxplots (stage, sex, tissue)")
message("    07-09: Dot plots with error bars (stage, sex, tissue)")
message("  DE PLOTS:")
message("    10: Volcano plots (one per comparison)")
message("    11: DE dot plot (family × comparison) - UPDATED colors")
message("  TABLE VISUALIZATIONS:")
message("    12: Summary heatmap (family × comparison)")
message("    13-14: Lollipop charts (DE counts)")
message("    15-16: Stacked bar charts (up/down direction)")
message("    17: UpSet plot (overlap between comparisons)")
message("  SAMPLE PLOTS:")
message("    18: PCA all samples - FIXED NA legend")
message("    18b: PCA all samples with cluster ellipses - NEW")
message("    19: PCA adults by tissue")
message("    19b: PCA adults by tissue with ellipses - NEW")
message("  GENE FAMILY HEATMAPS:")
message("    20: Individual family heatmaps (OBP, OR, CSP, P450, GST, UGT, SerPro, CAZyme)")
message("    20: Combined group heatmaps (Chemosensory, Detoxification, Digestion)")
message("\nData files:")
message("    visualization_data.RData")
message("    DE_family_summary_table.csv")

message("\n", paste(rep("=", 70), collapse = ""))
message("Pipeline completed: ", Sys.time())
message(paste(rep("=", 70), collapse = ""))
