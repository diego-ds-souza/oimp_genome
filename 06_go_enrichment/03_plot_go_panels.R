#!/usr/bin/env Rscript
#
# Figure 3: GO enrichment dot plots, one panel per comparative baseline and
# namespace. Point position is fold enrichment, colour is -log10(FDR) and size
# is the number of genes carrying the term.
#
# Terms are labelled by GO ID rather than by name, because the full names do
# not fit at the panel width used in the paper; the names are listed in the
# accompanying table.
#
# Run from the repository root:
#   Rscript 06_go_enrichment/03_plot_go_panels.R

# ------------------------------- settings ------------------------------------
RESULTS_DIR <- Sys.getenv("RESULTS_DIR", "results/06_go_enrichment/results")
file_suffix <- ".tsv"          # which enrichment tables to plot
TOP_N       <- 20              # top terms shown per namespace
MAX_FDR     <- 0.05            # FDR threshold
MIN_STUDY   <- 3               # minimum genes in the study set for a term
MIN_POP     <- 10              # minimum genes in the background for a term
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(scales)
})

if (!dir.exists(RESULTS_DIR)) stop("Results directory not found:\n  ", RESULTS_DIR)
setwd(RESULTS_DIR)

contrast_info <- list(
  expanded_vs_Aglabripennis = list(
    file   = paste0("GO_expanded_vs_Aglabripennis", file_suffix),
    label  = expression(italic("O. impluviata") * " vs " * italic("A. glabripennis")),
    prefix = "Fig3a"
  ),
  expanded_vs_Cerambycidae = list(
    file   = paste0("GO_expanded_vs_Cerambycidae", file_suffix),
    label  = expression(italic("O. impluviata") * " vs other Cerambycidae"),
    prefix = "Fig3b"
  ),
  expanded_vs_Chrysomelidae = list(
    file   = paste0("GO_expanded_vs_Chrysomelidae", file_suffix),
    label  = expression(italic("O. impluviata") * " vs Chrysomelidae"),
    prefix = "Fig3c"
  ),
  expanded_vs_Curculionoidea = list(
    file   = paste0("GO_expanded_vs_Curculionoidea", file_suffix),
    label  = expression(italic("O. impluviata") * " vs Curculionoidea"),
    prefix = "Fig3d"
  )
)

# Output panel names and order
panel_map <- c(
  "Cellular Component" = "CC",
  "Molecular Function" = "MF",
  "Biological process" = "BP"
)

# 1) Helpers for labelling and ordering terms.
std_cols <- function(nms) {
  nms <- trimws(nms)
  out <- character(length(nms))
  for (i in seq_along(nms)) {
    nm <- nms[i]
    if (nm == "GO")                       out[i] <- "GO"
    else if (grepl("^Term", nm))          out[i] <- "Term"
    else if (nm == "Name")                out[i] <- "Name"
    else if (nm == "Ratio in study")      out[i] <- "Ratio_study"
    else if (nm == "Ratio in pop")        out[i] <- "Ratio_pop"
    else if (nm == "fold_enrichment")     out[i] <- "fold_enrichment"
    else if (nm == "study_count")         out[i] <- "study_count"
    else if (nm == "pop_count")           out[i] <- "pop_count"
    else if (nm == "pvalue")              out[i] <- "pvalue"
    else if (nm == "p.adjust")            out[i] <- "p.adjust"
    else out[i] <- nm
  }
  out
}

split_ratio <- function(x) {
  s <- str_split_fixed(x, "/", 2)
  data.frame(a = as.numeric(s[,1]), b = as.numeric(s[,2]))
}

read_go_table <- function(fname) {
  if (!file.exists(fname)) {
    warning("File not found: ", fname)
    return(NULL)
  }
  df <- read_tsv(fname, col_types = cols(.default = col_character()), progress = FALSE)
  names(df) <- std_cols(names(df))
  
  needed <- c("GO", "Term", "Name", "pvalue", "p.adjust")
  if (!all(needed %in% names(df))) {
    warning("Skipping ", fname, " - missing columns: ", 
            paste(setdiff(needed, names(df)), collapse = ", "))
    return(NULL)
  }
  df
}

plot_one_namespace <- function(df_ns, title_expr, out_file) {
  # Use GO IDs instead of term names for better space utilization
  df_ns <- df_ns %>%
    mutate(
      TermLabel = GO,
      TermLabel = fct_reorder(TermLabel, RichFactor, .desc = FALSE)
    )
  
  # Breaks for GeneNumber legend
  rng <- range(df_ns$study_a, na.rm = TRUE)
  size_breaks <- pretty(rng, n = 4)
  size_breaks <- unique(size_breaks[size_breaks > 0])
  
  # Create plot
  p <- ggplot(df_ns, aes(x = RichFactor, y = TermLabel)) +
    geom_point(aes(size = study_a, color = neglog10Padj)) +
    scale_size_continuous(
      name = "Gene count",
      range = c(2, 8),
      breaks = size_breaks
    ) +
    scale_color_gradient(
      name = expression(-log[10] ~ "(FDR)"),
      low = "steelblue2",
      high = "firebrick2"
    ) +
    guides(
      color = guide_colorbar(order = 1),
      size = guide_legend(order = 2)
    ) +
    labs(
      x = "Fold enrichment",
      y = "GO terms",
      title = title_expr
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_line(linewidth = 0.3, colour = "grey90"),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      axis.text.y = element_text(size = 8),
      plot.title = element_text(hjust = 0, face = "plain", size = 10),
      plot.margin = margin(5, 10, 5, 5)
    )
  
  # Height scales with number of terms
  n_terms <- nrow(df_ns)
  plot_height <- max(4, 0.25 * n_terms + 1.5)
  
  ggsave(out_file, p, width = 6.5, height = plot_height)
  message("  Wrote: ", out_file, " (", n_terms, " terms)")
}

# Main loop.
message("\n", paste(rep("=", 70), collapse = ""))
message("GO Enrichment Visualization")
message(paste(rep("=", 70), collapse = ""), "\n")

for (key in names(contrast_info)) {
  info <- contrast_info[[key]]
  message("Processing: ", key)
  message("  File: ", info$file)
  
  df <- read_go_table(info$file)
  if (is.null(df)) next
  
  # Check if we have fold_enrichment column, otherwise calculate
  if ("fold_enrichment" %in% names(df)) {
    df <- df %>%
      mutate(
        study_a = as.numeric(study_count),
        pop_a = as.numeric(pop_count),
        RichFactor = as.numeric(fold_enrichment)
      )
  } else if ("study_count" %in% names(df)) {
    df <- df %>%
      mutate(
        study_a = as.numeric(study_count),
        study_total = as.numeric(study_total),
        pop_a = as.numeric(pop_count),
        pop_total = as.numeric(pop_total),
        RichFactor = ifelse(pop_a > 0, (study_a/study_total) / (pop_a/pop_total), NA_real_)
      )
  } else if ("Ratio_study" %in% names(df)) {
    rs <- split_ratio(df$Ratio_study)
    rp <- split_ratio(df$Ratio_pop)
    df <- df %>%
      mutate(
        study_a = rs$a,
        study_total = rs$b,
        pop_a = rp$a,
        pop_total = rp$b,
        RichFactor = ifelse(pop_a > 0, (study_a/study_total) / (pop_a/pop_total), NA_real_)
      )
  }
  
  # Calculate metrics
  df <- df %>%
    mutate(
      pvalue = as.numeric(pvalue),
      p.adjust = as.numeric(p.adjust),
      neglog10P = -log10(pvalue),
      neglog10Padj = -log10(p.adjust),
      Namespace = case_when(
        Term == "CC" ~ "Cellular Component",
        Term == "MF" ~ "Molecular Function",
        Term == "BP" ~ "Biological process",
        TRUE ~ Term
      )
    ) %>%
    filter(
      !is.na(p.adjust),
      p.adjust < MAX_FDR,
      study_a >= MIN_STUDY,
      pop_a >= MIN_POP,
      !is.na(RichFactor),
      !is.na(neglog10Padj)
    )
  
  if (nrow(df) == 0) {
    warning("No GO terms left after filtering for ", key, ". Skipping.")
    next
  }
  
  message("  Significant terms after filtering: ", nrow(df))
  
  # Process each namespace
  for (ns in names(panel_map)) {
    df_ns <- df %>%
      filter(Namespace == ns) %>%
      arrange(p.adjust) %>%
      slice_head(n = TOP_N)
    
    if (nrow(df_ns) == 0) {
      message("  [SKIP] No terms for namespace: ", ns)
      next
    }
    
    # Save term list used in that panel
    list_file <- paste0("list_GO_terms_", key, "_", panel_map[[ns]], ".tsv")
    write_tsv(df_ns, list_file)
    
    # Build title expression
    title_expr <- bquote(.(info$label) ~ " | " ~ .(ns))
    
    # Output file
    out_file <- paste0(info$prefix, "_", panel_map[[ns]], "_GO_enrichment_", key, ".pdf")
    plot_one_namespace(df_ns, title_expr, out_file)
  }
  
  message("")
}

message(paste(rep("=", 70), collapse = ""))
message("Done!")
message(paste(rep("=", 70), collapse = ""))
