#!/usr/bin/env Rscript
#
# Chromosome-scale organization of the host-interface gene families (Fig. 7).
#
# Gene coordinates are read from the GFF3 and converted from scaffold-relative
# to chromosome-scale positions using the scaffold-to-chromosome assignments
# from the synteny analysis. Scaffold orientation comes from the telomere
# calls: a scaffold whose telomeres imply the opposite polarity is reversed
# before its coordinates are made cumulative, so every chromosome has a
# consistent direction.
#
# From those positions the script derives gene density in sliding windows,
# tandem clusters, and superloci, then draws the ideograms.
#
# Run from the repository root:
#   Rscript 08_chromosomal_organization/01_chromosome_ideograms.R

# ------------------------------- settings ------------------------------------
GFF3_FILE         <- Sys.getenv("GFF3_FILE",
  "results/03_annotation/curated/Oncideres_impluviata.clean_AGAT.renamed.products.fixed.renumbered.gff3")
SCAFFOLD_MAP_FILE <- "08_chromosomal_organization/scaffolds_to_chromosome.csv"
GENE_NAMES_FILE   <- Sys.getenv("GENE_NAMES_FILE",
                                "results/07_transcriptomics/ref/gene_names.tsv")
FAMILIES_DIR      <- Sys.getenv("FAMILIES_DIR", "results/07_transcriptomics/families")
DESEQ_RESULTS     <- Sys.getenv("DESEQ_RESULTS",
                                "results/07_transcriptomics/deseq2/deseq2_results.RData")
OUTPUT_DIR        <- Sys.getenv("OUTPUT_DIR", "results/08_chromosomal_organization")
PLOTS_DIR         <- file.path(OUTPUT_DIR, "plots")

MAX_CLUSTER_DISTANCE <- 100000   # bp between consecutive genes of one cluster
MIN_CLUSTER_GENES    <- 2        # genes required to call a cluster
SUPERLOCUS_DISTANCE  <- 50000    # bp between clusters of a superlocus
DENSITY_WINDOW       <- 500000   # sliding window for gene density
# -----------------------------------------------------------------------------

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

# 1) Load packages.
message("\n", paste(rep("=", 70), collapse = ""))
message("Chromosome Ideogram with Synteny-Based Mapping (UPDATED)")
message(paste(rep("=", 70), collapse = ""))
message("\nLoading packages...")

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (has_patchwork) library(patchwork)

# 2) Load scaffold-to-chromosome mapping with orientation.
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 1: Loading Scaffold-to-Chromosome Mapping")
message(paste(rep("=", 70), collapse = ""))

message("\n[Step 1] Loading synteny-based chromosome mapping...")

if (!file.exists(SCAFFOLD_MAP_FILE)) {
  stop("Scaffold mapping file not found: ", SCAFFOLD_MAP_FILE)
}

scaffold_map <- read_csv(SCAFFOLD_MAP_FILE, show_col_types = FALSE)

# Clean up and determine scaffold orientation
scaffold_map <- scaffold_map %>%
  mutate(
    chromosome = ifelse(chromosome == "" | is.na(chromosome), NA_character_, chromosome),
    has_chromosome = !is.na(chromosome),
    telomere_5prime = tolower(telomere_5prime) == "yes",
    telomere_3prime = tolower(telomere_3prime) == "yes"
  )

message("  Scaffolds with chromosome assignment: ", sum(scaffold_map$has_chromosome))

# -----------------------------------------------------------------------------
# Determine proper scaffold order and orientation within each chromosome
# Based on telomere positions:
# - First scaffold should have 5' telomere at chromosome start
# - Last scaffold should have 3' telomere at chromosome end
# - If scaffold has wrong telomere, it needs to be reverse-complemented
# -----------------------------------------------------------------------------

message("\n[Step 2] Determining scaffold orientations...")

# For chromosomes with multiple scaffolds, determine order and orientation
scaffold_order <- scaffold_map %>%
  filter(has_chromosome) %>%
  group_by(chromosome) %>%
  mutate(n_scaffolds = n()) %>%
  ungroup()

# Manually define scaffold order and orientation based on telomere analysis
# From your data:
# chr_1: scaffold_9 (5' tel) + scaffold_5 (3' tel) - CORRECT ORDER
# chr_2: scaffold_2 (5' tel) + scaffold_12 (5' tel) - scaffold_12 needs REVERSAL
# chr_6: scaffold_10 (5' tel) + scaffold_13 (3' tel) - CORRECT ORDER

scaffold_orientation <- tribble(
  ~scaffold, ~chromosome, ~order_in_chr, ~is_reversed,
  # chr_1: scaffold_9 first (has 5' tel), scaffold_5 second (has 3' tel)
  "scaffold_9", "chr_1", 1, FALSE,
  "scaffold_5", "chr_1", 2, FALSE,
  # chr_2: scaffold_2 first (has 5' tel), scaffold_12 second (has 5' tel, so REVERSE it)
  "scaffold_2", "chr_2", 1, FALSE,
  "scaffold_12", "chr_2", 2, TRUE,  # Reversed so 5' becomes 3' at chr end
  # chr_3: single scaffold with both telomeres
  "scaffold_1", "chr_3", 1, FALSE,
  # chr_4: single scaffold, no telomeres
  "scaffold_3", "chr_4", 1, FALSE,
  # chr_5: single scaffold with both telomeres
  "scaffold_4", "chr_5", 1, FALSE,
  # chr_6: scaffold_10 first (has 5' tel), scaffold_13 second (has 3' tel)
  "scaffold_10", "chr_6", 1, FALSE,
  "scaffold_13", "chr_6", 2, FALSE,
  # chr_7: single scaffold (has 5' tel only)
  "scaffold_6", "chr_7", 1, FALSE,
  # chr_8: single scaffold with both telomeres
  "scaffold_7", "chr_8", 1, FALSE,
  # chr_9: single scaffold with both telomeres
  "scaffold_8", "chr_9", 1, FALSE,
  # chr_10: single scaffold with both telomeres
  "scaffold_11", "chr_10", 1, FALSE
)

# Join with scaffold lengths
scaffold_orientation <- scaffold_orientation %>%
  left_join(scaffold_map %>% select(scaffold, length_bp, telomere_5prime, telomere_3prime), 
            by = "scaffold")

# Calculate cumulative offsets for each scaffold within chromosome
scaffold_orientation <- scaffold_orientation %>%
  arrange(chromosome, order_in_chr) %>%
  group_by(chromosome) %>%
  mutate(
    offset = cumsum(lag(length_bp, default = 0)),
    chr_start = offset,
    chr_end = offset + length_bp
  ) %>%
  ungroup()

message("\n  Scaffold orientations:")
print(scaffold_orientation %>% select(chromosome, scaffold, order_in_chr, is_reversed, length_bp, offset))

# Calculate chromosome lengths
chr_lengths <- scaffold_orientation %>%
  group_by(chromosome) %>%
  summarise(
    length = sum(length_bp),
    n_scaffolds = n(),
    scaffolds = paste(scaffold, collapse = " + "),
    .groups = "drop"
  ) %>%
  mutate(chr_num = as.numeric(str_extract(chromosome, "\\d+"))) %>%
  arrange(chr_num)

message("\n  Chromosome lengths:")
print(chr_lengths)

# 3) Parse GFF3 and map genes to chromosomes.
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 2: Parsing GFF3 and Mapping Genes")
message(paste(rep("=", 70), collapse = ""))

# Function to parse GFF3
parse_gff3_genes <- function(gff3_file) {
  message("\n[Step 3] Reading GFF3 file: ", basename(gff3_file))
  
  gff_lines <- readLines(gff3_file)
  gff_lines <- gff_lines[!grepl("^#", gff_lines)]
  gff_lines <- gff_lines[gff_lines != ""]
  
  gff_df <- read.delim(
    text = gff_lines, 
    header = FALSE, 
    sep = "\t",
    stringsAsFactors = FALSE,
    col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  )
  
  genes_df <- gff_df %>%
    filter(type == "gene") %>%
    mutate(
      gene_id = str_extract(attributes, "ID=[^;]+") %>% str_remove("ID="),
      gene_name = str_extract(attributes, "Name=[^;]+") %>% str_remove("Name=")
    ) %>%
    select(seqid, start, end, strand, gene_id, gene_name) %>%
    rename(scaffold = seqid) %>%
    arrange(scaffold, start)
  
  message("  Genes found: ", nrow(genes_df))
  return(genes_df)
}

# Load gene positions
gene_positions_file <- file.path(OUTPUT_DIR, "gene_positions.tsv")

if (file.exists(GFF3_FILE)) {
  gene_positions <- parse_gff3_genes(GFF3_FILE)
  write_tsv(gene_positions, gene_positions_file)
} else if (file.exists(gene_positions_file)) {
  message("\n[Step 3] Loading pre-extracted gene positions...")
  gene_positions <- read_tsv(gene_positions_file, show_col_types = FALSE)
  message("  Genes loaded: ", nrow(gene_positions))
} else {
  stop("GFF3 file not found and no pre-extracted positions available.")
}

# 4) Map genes to chromosome coordinates (handling reversed scaffolds).
message("\n[Step 4] Mapping genes to chromosome coordinates...")

gene_positions_chr <- gene_positions %>%
  left_join(
    scaffold_orientation %>% select(scaffold, chromosome, offset, length_bp, is_reversed),
    by = "scaffold"
  ) %>%
  filter(!is.na(chromosome)) %>%
  mutate(
    # For reversed scaffolds, flip coordinates
    chr_start = ifelse(is_reversed,
                       offset + (length_bp - end),
                       offset + start),
    chr_end = ifelse(is_reversed,
                     offset + (length_bp - start),
                     offset + end),
    # Also flip strand for reversed scaffolds
    chr_strand = ifelse(is_reversed,
                        ifelse(strand == "+", "-", "+"),
                        strand)
  )

message("  Genes mapped to chromosomes: ", nrow(gene_positions_chr))

# 5) Calculate gene density.
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 3: Calculating Gene Density")
message(paste(rep("=", 70), collapse = ""))

message("\n[Step 5] Calculating gene density in ", DENSITY_WINDOW/1000, " kb windows...")

# Create windows for each chromosome
density_list <- list()
for (i in 1:nrow(chr_lengths)) {
  chr <- chr_lengths[i, ]
  windows <- data.frame(
    chromosome = chr$chromosome,
    chr_num = chr$chr_num,
    window_start = seq(0, chr$length - 1, by = DENSITY_WINDOW),
    stringsAsFactors = FALSE
  )
  windows$window_end <- pmin(windows$window_start + DENSITY_WINDOW, chr$length)
  windows$window_mid <- (windows$window_start + windows$window_end) / 2
  density_list[[i]] <- windows
}
density_data <- bind_rows(density_list)

# Count genes per window
density_data$n_genes <- sapply(1:nrow(density_data), function(i) {
  sum(gene_positions_chr$chromosome == density_data$chromosome[i] & 
      gene_positions_chr$chr_start >= density_data$window_start[i] & 
      gene_positions_chr$chr_start < density_data$window_end[i])
})

density_data$density <- density_data$n_genes / (DENSITY_WINDOW / 1e6)  # genes per Mb

message("  Gene density calculated for ", nrow(density_data), " windows")
message("  Mean density: ", round(mean(density_data$density), 1), " genes/Mb")
message("  Range: ", round(min(density_data$density), 1), " - ", 
        round(max(density_data$density), 1), " genes/Mb")

# 6) Load gene family annotations.
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 4: Loading Gene Family Annotations")
message(paste(rep("=", 70), collapse = ""))

message("\n[Step 6] Loading gene annotations...")

if (file.exists(GENE_NAMES_FILE)) {
  gene_names <- read_tsv(GENE_NAMES_FILE, show_col_types = FALSE)
  message("  Gene names loaded: ", nrow(gene_names))
} else {
  gene_names <- NULL
}

# Define and load gene families
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

# Join with chromosome positions
gene_positions_annotated <- gene_positions_chr %>%
  left_join(family_assignments, by = "gene_id") %>%
  left_join(gene_names %>% select(gene_id, name, product), by = "gene_id")

write_tsv(gene_positions_annotated, file.path(OUTPUT_DIR, "gene_positions_chromosome_level.tsv"))

# 7) Identify gene clusters.
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 5: Identifying Gene Clusters (Tandem Arrays)")
message(paste(rep("=", 70), collapse = ""))

message("\n[Step 7] Identifying gene clusters...")
message("  Definition: ≥", MIN_CLUSTER_GENES, " genes from same family within ", 
        MAX_CLUSTER_DISTANCE/1000, " kb")

identify_clusters <- function(positions_df, family_col = "family", 
                              max_distance = MAX_CLUSTER_DISTANCE, 
                              min_genes = MIN_CLUSTER_GENES) {
  clusters_list <- list()
  cluster_id <- 0
  
  families <- unique(positions_df[[family_col]])
  families <- families[!is.na(families)]
  
  for (fam in families) {
    fam_genes <- positions_df %>%
      filter(!!sym(family_col) == fam & !is.na(chromosome)) %>%
      arrange(chromosome, chr_start)
    
    if (nrow(fam_genes) < min_genes) next
    
    for (chr in unique(fam_genes$chromosome)) {
      chr_genes <- fam_genes %>%
        filter(chromosome == chr) %>%
        arrange(chr_start)
      
      if (nrow(chr_genes) < min_genes) next
      
      current_cluster <- c(1)
      
      for (i in 2:nrow(chr_genes)) {
        dist_to_prev <- chr_genes$chr_start[i] - chr_genes$chr_end[i-1]
        
        if (dist_to_prev <= max_distance) {
          current_cluster <- c(current_cluster, i)
        } else {
          if (length(current_cluster) >= min_genes) {
            cluster_id <- cluster_id + 1
            cluster_genes <- chr_genes[current_cluster, ]
            
            clusters_list[[cluster_id]] <- data.frame(
              cluster_id = cluster_id,
              family = fam,
              group = unique(cluster_genes$group),
              chromosome = chr,
              scaffold = paste(unique(cluster_genes$scaffold), collapse = ","),
              start = min(cluster_genes$chr_start),
              end = max(cluster_genes$chr_end),
              n_genes = length(current_cluster),
              genes = paste(cluster_genes$gene_id, collapse = ","),
              span_kb = round((max(cluster_genes$chr_end) - min(cluster_genes$chr_start)) / 1000, 1)
            )
          }
          current_cluster <- c(i)
        }
      }
      
      # Last cluster
      if (length(current_cluster) >= min_genes) {
        cluster_id <- cluster_id + 1
        cluster_genes <- chr_genes[current_cluster, ]
        
        clusters_list[[cluster_id]] <- data.frame(
          cluster_id = cluster_id,
          family = fam,
          group = unique(cluster_genes$group),
          chromosome = chr,
          scaffold = paste(unique(cluster_genes$scaffold), collapse = ","),
          start = min(cluster_genes$chr_start),
          end = max(cluster_genes$chr_end),
          n_genes = length(current_cluster),
          genes = paste(cluster_genes$gene_id, collapse = ","),
          span_kb = round((max(cluster_genes$chr_end) - min(cluster_genes$chr_start)) / 1000, 1)
        )
      }
    }
  }
  
  if (length(clusters_list) > 0) {
    return(bind_rows(clusters_list))
  } else {
    return(data.frame())
  }
}

gene_clusters <- identify_clusters(gene_positions_annotated)

if (nrow(gene_clusters) > 0) {
  message("  Total clusters found: ", nrow(gene_clusters))
  message("  Total clustered genes: ", sum(gene_clusters$n_genes))
  
  write_tsv(gene_clusters, file.path(OUTPUT_DIR, "gene_clusters_chromosome_level.tsv"))
  
  message("\n  Top 10 largest clusters:")
  print(gene_clusters %>% arrange(desc(n_genes)) %>% head(10) %>% 
          select(family, chromosome, n_genes, span_kb))
}

# 8) Identify superloci (co-localized clusters from different families).
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 6: Identifying Superloci")
message(paste(rep("=", 70), collapse = ""))

message("\n[Step 8] Identifying superloci (co-localized clusters)...")

# Find clusters that overlap or are very close
identify_superloci <- function(clusters_df, max_distance = SUPERLOCUS_DISTANCE) {
  superloci <- list()
  superlocus_id <- 0
  
  for (chr in unique(clusters_df$chromosome)) {
    chr_clusters <- clusters_df %>%
      filter(chromosome == chr) %>%
      arrange(start)
    
    if (nrow(chr_clusters) < 2) next
    
    i <- 1
    while (i < nrow(chr_clusters)) {
      current_superlocus <- chr_clusters[i, ]
      current_end <- current_superlocus$end
      
      j <- i + 1
      while (j <= nrow(chr_clusters)) {
        next_cluster <- chr_clusters[j, ]
        
        # Check if clusters overlap or are within max_distance
        if (next_cluster$start <= current_end + max_distance) {
          # Different family? This is a superlocus!
          if (next_cluster$family != current_superlocus$family[1]) {
            current_superlocus <- rbind(current_superlocus, next_cluster)
            current_end <- max(current_end, next_cluster$end)
          }
          j <- j + 1
        } else {
          break
        }
      }
      
      # If we found multiple families co-localized
      if (length(unique(current_superlocus$family)) > 1) {
        superlocus_id <- superlocus_id + 1
        superloci[[superlocus_id]] <- data.frame(
          superlocus_id = superlocus_id,
          chromosome = chr,
          start = min(current_superlocus$start),
          end = max(current_superlocus$end),
          n_clusters = nrow(current_superlocus),
          families = paste(unique(current_superlocus$family), collapse = " + "),
          total_genes = sum(current_superlocus$n_genes),
          span_kb = round((max(current_superlocus$end) - min(current_superlocus$start)) / 1000, 1)
        )
      }
      
      i <- j
    }
  }
  
  if (length(superloci) > 0) {
    return(bind_rows(superloci))
  } else {
    return(data.frame())
  }
}

superloci <- identify_superloci(gene_clusters)

if (nrow(superloci) > 0) {
  message("  Superloci found: ", nrow(superloci))
  print(superloci)
  write_tsv(superloci, file.path(OUTPUT_DIR, "superloci.tsv"))
  
  # Highlight the chr_3 digestion superlocus
  chr3_superlocus <- superloci %>%
    filter(chromosome == "chr_3" & grepl("UGT|CAZyme", families))
  
  if (nrow(chr3_superlocus) > 0) {
    message("\n  *** Chr 3 Digestion Superlocus identified ***")
    print(chr3_superlocus)
  }
} else {
  message("  No superloci found")
  superloci <- data.frame()
}

# 9) Load DE results.
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 7: Loading DE Results")
message(paste(rep("=", 70), collapse = ""))

if (file.exists(DESEQ_RESULTS)) {
  load(DESEQ_RESULTS)
  message("  Comparisons loaded: ", length(all_results))
  
  comparison_labels <- c(
    "antennae_MvF" = "M vs F (Ant)",
    "thorax_MvF" = "M vs F (Th)",
    "adults_MvF" = "M vs F (All)",
    "tissue" = "Ant vs Th",
    "larva_vs_male" = "L vs M",
    "larva_vs_female" = "L vs F",
    "larva_vs_adults" = "L vs A"
  )
  
  top_de_genes_list <- list()
  
  for (comp_name in names(all_results)) {
    de_res <- as.data.frame(all_results[[comp_name]]$results)
    de_res$gene_id <- rownames(de_res)
    
    top_genes <- de_res %>%
      filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1) %>%
      arrange(padj) %>%
      head(10) %>%
      mutate(comparison = comp_name)
    
    if (nrow(top_genes) > 0) {
      top_de_genes_list[[comp_name]] <- top_genes
    }
  }
  
  top_de_genes <- bind_rows(top_de_genes_list) %>%
    left_join(gene_positions_chr %>% select(gene_id, chromosome, chr_start, chr_end), by = "gene_id") %>%
    left_join(gene_names %>% select(gene_id, name, product), by = "gene_id")
  
  write_tsv(top_de_genes, file.path(OUTPUT_DIR, "top_DE_genes_chromosome_level.tsv"))
} else {
  message("  DESeq2 results not found")
  all_results <- NULL
  top_de_genes <- NULL
}

# 10) Create chromosome ideogram figures.
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 8: Creating Chromosome Ideogram Figures")
message(paste(rep("=", 70), collapse = ""))

# 11) Define colors.
family_colors <- c(
  "OBP" = "#1B9E77", "CSP" = "#66C2A5", "OR" = "#B2E2E2",
  "GR" = "#238B45", "IR" = "#74C476", "SNMP" = "#BAE4B3",
  "P450" = "#D95F02", "GST" = "#FC8D62", "UGT" = "#FDBF6F",
  "CCE" = "#E31A1C", "ABC" = "#FB9A99",
  "SerPro" = "#1F78B4", "CysPro" = "#6A3D9A", "AspPro" = "#A6CEE3",
  "MetPro" = "#CAB2D6", "CAZyme" = "#33A02C"
)

# Prepare data
chr_data <- chr_lengths %>% arrange(chr_num)

family_positions_chr <- gene_positions_annotated %>%
  filter(!is.na(chromosome) & !is.na(family)) %>%
  left_join(chr_lengths %>% select(chromosome, chr_num), by = "chromosome")

# Prepare scaffold junction data
scaffold_junctions <- scaffold_orientation %>%
  group_by(chromosome) %>%
  filter(n() > 1) %>%
  filter(order_in_chr > 1) %>%
  ungroup() %>%
  left_join(chr_lengths %>% select(chromosome, chr_num), by = "chromosome")

# 12) Helper function to create scaffold junction markers (double diagonal lines).
create_junction_geom <- function(junctions_df, y_offset = 0.35, gap_width = 0.02) {
  if (nrow(junctions_df) == 0) return(NULL)
  
  # Create data for diagonal lines
  junction_lines <- junctions_df %>%
    rowwise() %>%
    do({
      j <- .
      x_center <- j$offset
      y_center <- j$chr_num
      
      # Create two parallel diagonal lines
      data.frame(
        chromosome = j$chromosome,
        chr_num = j$chr_num,
        line_id = c(1, 1, 2, 2),
        x = c(x_center - gap_width * j$offset * 0.001 - 500000, 
              x_center - gap_width * j$offset * 0.001 + 500000,
              x_center + gap_width * j$offset * 0.001 - 500000, 
              x_center + gap_width * j$offset * 0.001 + 500000),
        y = c(y_center - y_offset - 0.1, y_center + y_offset + 0.1,
              y_center - y_offset - 0.1, y_center + y_offset + 0.1)
      )
    }) %>%
    ungroup()
  
  junction_lines
}

# 13) Main ideogram with all families, scaffold junctions, and density.
message("\n[Step 9] Creating main ideogram with gene density...")

# Calculate y positions
chr_y_base <- 1:10
density_y_offset <- 0.6  # Height for density track

p_main <- ggplot() +
  # Gene density track (as area/ribbon above chromosomes)
  geom_ribbon(
    data = density_data,
    aes(x = window_mid, 
        ymin = chr_num + 0.4, 
        ymax = chr_num + 0.4 + (density / max(density_data$density)) * 0.5,
        group = chromosome),
    fill = "grey60", alpha = 0.5
  ) +
  # Chromosome backbones
  geom_rect(
    data = chr_data,
    aes(xmin = 0, xmax = length, 
        ymin = chr_num - 0.35, ymax = chr_num + 0.35),
    fill = "grey20", color = "grey10", linewidth = 0.3
  ) +
  # Scaffold junctions (white gap with diagonal lines)
  geom_rect(
    data = scaffold_junctions,
    aes(xmin = offset - 800000, xmax = offset + 800000,
        ymin = chr_num - 0.4, ymax = chr_num + 0.4),
    fill = "white", color = NA
  ) +
  # Junction diagonal lines
  geom_segment(
    data = scaffold_junctions,
    aes(x = offset - 600000, xend = offset - 200000,
        y = chr_num - 0.4, yend = chr_num + 0.4),
    color = "grey30", linewidth = 0.8
  ) +
  geom_segment(
    data = scaffold_junctions,
    aes(x = offset + 200000, xend = offset + 600000,
        y = chr_num - 0.4, yend = chr_num + 0.4),
    color = "grey30", linewidth = 0.8
  ) +
  # Telomere markers
  geom_point(
    data = scaffold_orientation %>% 
      filter(telomere_5prime & order_in_chr == 1) %>%
      left_join(chr_lengths %>% select(chromosome, chr_num), by = "chromosome"),
    aes(x = 0, y = chr_num),
    color = "#E41A1C", size = 3
  ) +
  geom_point(
    data = scaffold_orientation %>% 
      group_by(chromosome) %>%
      filter(order_in_chr == max(order_in_chr)) %>%
      filter(telomere_3prime | (is_reversed & telomere_5prime)) %>%
      ungroup() %>%
      left_join(chr_lengths %>% select(chromosome, chr_num, length), by = "chromosome"),
    aes(x = length, y = chr_num),
    color = "#E41A1C", size = 3
  ) +
  # Gene family markers
  geom_segment(
    data = family_positions_chr,
    aes(x = chr_start, xend = chr_start,
        y = chr_num - 0.28, yend = chr_num + 0.28,
        color = family),
    linewidth = 0.3, alpha = 0.7
  ) +
  scale_color_manual(values = family_colors, name = "Gene Family") +
  scale_x_continuous(
    labels = label_number(scale = 1e-6, suffix = " Mb"),
    expand = c(0.02, 0)
  ) +
  scale_y_continuous(
    breaks = chr_data$chr_num,
    labels = paste0("chr ", chr_data$chr_num),
    expand = c(0.08, 0)
  ) +
  labs(
    title = expression(paste("Chromosomal Distribution of Gene Families in ", italic("Oncideres impluviata"))),
    subtitle = "Based on synteny with M. alternatus | Red dots = telomeres | Grey ribbon = gene density | Double lines = scaffold junctions",
    x = "Position",
    y = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    axis.text.y = element_text(size = 11, face = "bold"),
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(ncol = 1, override.aes = list(linewidth = 2, alpha = 1)))

ggsave(file.path(PLOTS_DIR, "05_chromosome_ideogram_synteny_all_families.pdf"), 
       p_main, width = 16, height = 11)
message("  Saved: 05_chromosome_ideogram_synteny_all_families.pdf")

# 14) Ideograms by functional group.
message("\n[Step 10] Creating ideograms by functional group...")

for (group_name in names(FAMILY_GROUPS)) {
  group_families <- FAMILY_GROUPS[[group_name]]
  group_positions <- family_positions_chr %>%
    filter(family %in% group_families)
  
  if (nrow(group_positions) == 0) next
  
  p_group <- ggplot() +
    # Gene density track
    geom_ribbon(
      data = density_data,
      aes(x = window_mid, 
          ymin = chr_num + 0.4, 
          ymax = chr_num + 0.4 + (density / max(density_data$density)) * 0.5,
          group = chromosome),
      fill = "grey70", alpha = 0.4
    ) +
    # Chromosome backbones
    geom_rect(
      data = chr_data,
      aes(xmin = 0, xmax = length, 
          ymin = chr_num - 0.35, ymax = chr_num + 0.35),
      fill = "grey85", color = "grey60", linewidth = 0.3
    ) +
    # Scaffold junctions
    geom_rect(
      data = scaffold_junctions,
      aes(xmin = offset - 800000, xmax = offset + 800000,
          ymin = chr_num - 0.4, ymax = chr_num + 0.4),
      fill = "white", color = NA
    ) +
    geom_segment(
      data = scaffold_junctions,
      aes(x = offset - 600000, xend = offset - 200000,
          y = chr_num - 0.4, yend = chr_num + 0.4),
      color = "grey50", linewidth = 0.8
    ) +
    geom_segment(
      data = scaffold_junctions,
      aes(x = offset + 200000, xend = offset + 600000,
          y = chr_num - 0.4, yend = chr_num + 0.4),
      color = "grey50", linewidth = 0.8
    ) +
    # Telomeres
    geom_point(
      data = scaffold_orientation %>% 
        filter(telomere_5prime & order_in_chr == 1) %>%
        left_join(chr_lengths %>% select(chromosome, chr_num), by = "chromosome"),
      aes(x = 0, y = chr_num),
      color = "#E41A1C", size = 3
    ) +
    geom_point(
      data = scaffold_orientation %>% 
        group_by(chromosome) %>%
        filter(order_in_chr == max(order_in_chr)) %>%
        filter(telomere_3prime | (is_reversed & telomere_5prime)) %>%
        ungroup() %>%
        left_join(chr_lengths %>% select(chromosome, chr_num, length), by = "chromosome"),
      aes(x = length, y = chr_num),
      color = "#E41A1C", size = 3
    ) +
    # Gene markers
    geom_segment(
      data = group_positions,
      aes(x = chr_start, xend = chr_start,
          y = chr_num - 0.28, yend = chr_num + 0.28,
          color = family),
      linewidth = 0.5, alpha = 0.8
    ) +
    scale_color_manual(values = family_colors[group_families], name = "Gene Family") +
    scale_x_continuous(
      labels = label_number(scale = 1e-6, suffix = " Mb"),
      expand = c(0.02, 0)
    ) +
    scale_y_continuous(
      breaks = chr_data$chr_num,
      labels = paste0("chr ", chr_data$chr_num),
      expand = c(0.08, 0)
    ) +
    labs(
      title = paste0(group_name, " Genes — Chromosomal Distribution"),
      subtitle = paste0(nrow(group_positions), " genes | Grey ribbon = gene density | Double lines = scaffold junctions"),
      x = "Position",
      y = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      axis.text.y = element_text(size = 11, face = "bold"),
      legend.position = "right",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Add superlocus highlight for digestion group
  if (group_name == "Digestion" && nrow(superloci) > 0) {
    digestion_superloci <- superloci %>%
      filter(grepl("UGT|CAZyme|SerPro", families)) %>%
      left_join(chr_lengths %>% select(chromosome, chr_num), by = "chromosome")
    
    if (nrow(digestion_superloci) > 0) {
      p_group <- p_group +
        # Highlight superlocus with bracket
        geom_rect(
          data = digestion_superloci,
          aes(xmin = start - 500000, xmax = end + 500000,
              ymin = chr_num - 0.5, ymax = chr_num + 0.5),
          fill = NA, color = "#FF6B6B", linewidth = 1.5, linetype = "dashed"
        ) +
        # Label
        geom_label(
          data = digestion_superloci,
          aes(x = (start + end) / 2, y = chr_num + 1.0, 
              label = paste0("Superlocus\n", families)),
          size = 2.5, fill = "#FFE5E5", color = "#CC0000", fontface = "bold"
        )
    }
  }
  
  filename <- paste0("06_chromosome_ideogram_synteny_", tolower(group_name), ".pdf")
  ggsave(file.path(PLOTS_DIR, filename), p_group, width = 16, height = 11)
  message("  Saved: ", filename)
}

# 15) Cluster zoom plots with improved labeling.
if (nrow(gene_clusters) > 0) {
  message("\n[Step 11] Creating cluster zoom plots...")
  
  top_cluster_ids <- gene_clusters %>%
    arrange(desc(n_genes)) %>%
    head(6) %>%
    pull(cluster_id)
  
  cluster_plots <- list()
  
  for (i in seq_along(top_cluster_ids)) {
    clust <- gene_clusters %>% filter(cluster_id == top_cluster_ids[i])
    cluster_gene_ids <- unlist(strsplit(clust$genes, ","))
    
    cluster_genes <- gene_positions_annotated %>%
      filter(gene_id %in% cluster_gene_ids) %>%
      arrange(chr_start) %>%
      mutate(
        gene_num = row_number(),
        # Shorten gene ID for display
        short_id = str_replace(gene_id, "FUN_0*", "")
      )
    
    # Stagger labels to avoid overlap
    cluster_genes <- cluster_genes %>%
      mutate(
        label_y = ifelse(gene_num %% 3 == 1, 0.7, 
                         ifelse(gene_num %% 3 == 2, 0.5, 0.9))
      )
    
    p_zoom <- ggplot(cluster_genes) +
      # Gene rectangles with arrows showing strand
      geom_rect(
        aes(xmin = chr_start, xmax = chr_end, 
            ymin = -0.25, ymax = 0.25,
            fill = chr_strand),
        color = "black", linewidth = 0.3
      ) +
      # Gene labels (staggered)
      geom_text(
        aes(x = (chr_start + chr_end) / 2, y = label_y, label = short_id),
        angle = 60, hjust = 0, size = 2, color = "grey30"
      ) +
      # Connecting lines to labels
      geom_segment(
        aes(x = (chr_start + chr_end) / 2, xend = (chr_start + chr_end) / 2,
            y = 0.28, yend = label_y - 0.05),
        color = "grey70", linewidth = 0.3
      ) +
      scale_fill_manual(
        values = c("+" = family_colors[clust$family], "-" = "white"),
        labels = c("+" = "Sense (+)", "-" = "Antisense (−)"),
        name = "Strand"
      ) +
      scale_x_continuous(
        labels = label_number(scale = 1e-3, suffix = " kb")
      ) +
      coord_cartesian(ylim = c(-0.4, 1.3)) +
      labs(
        title = paste0(clust$family, " cluster — ", clust$chromosome, " (", clust$n_genes, " genes)"),
        subtitle = paste0("Span: ", clust$span_kb, " kb | Scaffold: ", clust$scaffold, 
                          " | Position: ", round(clust$start/1e6, 2), "-", round(clust$end/1e6, 2), " Mb"),
        x = "Chromosomal Position"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey50"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom"
      )
    
    cluster_plots[[i]] <- p_zoom
  }
  
  if (has_patchwork && length(cluster_plots) >= 2) {
    combined <- wrap_plots(cluster_plots, ncol = 2) +
      plot_annotation(
        title = "Gene Cluster Zoom Views — Largest Tandem Arrays",
        subtitle = "Colored = sense strand (+) | White = antisense strand (−) | Numbers = gene IDs (FUN_ prefix removed)",
        theme = theme(
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10, color = "grey40")
        )
      )
    ggsave(file.path(PLOTS_DIR, "07_gene_clusters_zoom_synteny.pdf"), 
           combined, width = 16, height = 16)
    message("  Saved: 07_gene_clusters_zoom_synteny.pdf")
  }
}

# 16) Top DE genes ideogram.
if (!is.null(top_de_genes) && nrow(top_de_genes) > 0) {
  message("\n[Step 12] Creating top DE genes ideogram...")
  
  top_de_positions <- top_de_genes %>%
    filter(!is.na(chromosome)) %>%
    left_join(chr_lengths %>% select(chromosome, chr_num), by = "chromosome")
  
  p_de <- ggplot() +
    geom_rect(
      data = chr_data,
      aes(xmin = 0, xmax = length, 
          ymin = chr_num - 0.35, ymax = chr_num + 0.35),
      fill = "grey90", color = "grey60", linewidth = 0.3
    ) +
    # Scaffold junctions
    geom_rect(
      data = scaffold_junctions,
      aes(xmin = offset - 800000, xmax = offset + 800000,
          ymin = chr_num - 0.4, ymax = chr_num + 0.4),
      fill = "white", color = NA
    ) +
    geom_segment(
      data = scaffold_junctions,
      aes(x = offset - 600000, xend = offset - 200000,
          y = chr_num - 0.4, yend = chr_num + 0.4),
      color = "grey50", linewidth = 0.8
    ) +
    geom_segment(
      data = scaffold_junctions,
      aes(x = offset + 200000, xend = offset + 600000,
          y = chr_num - 0.4, yend = chr_num + 0.4),
      color = "grey50", linewidth = 0.8
    ) +
    geom_point(
      data = top_de_positions,
      aes(x = chr_start, y = chr_num, color = comparison),
      size = 3, alpha = 0.8
    ) +
    scale_color_brewer(palette = "Set1", name = "Comparison",
                       labels = comparison_labels) +
    scale_x_continuous(
      labels = label_number(scale = 1e-6, suffix = " Mb"),
      expand = c(0.02, 0)
    ) +
    scale_y_continuous(
      breaks = chr_data$chr_num,
      labels = paste0("chr ", chr_data$chr_num),
      expand = c(0.08, 0)
    ) +
    labs(
      title = "Chromosomal Positions of Top Differentially Expressed Genes",
      subtitle = "Top 10 most significant genes from each comparison | Double lines = scaffold junctions",
      x = "Position",
      y = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 11, face = "bold"),
      legend.position = "right",
      panel.grid.major.y = element_blank()
    )
  
  ggsave(file.path(PLOTS_DIR, "08_chromosome_ideogram_top_DE_genes.pdf"), 
         p_de, width = 16, height = 10)
  message("  Saved: 08_chromosome_ideogram_top_DE_genes.pdf")
}

# 17) Summary tables.
message("\n", paste(rep("=", 70), collapse = ""))
message("PART 9: Creating Summary Tables")
message(paste(rep("=", 70), collapse = ""))

# Chromosome summary with clusters
chr_summary <- chr_lengths %>%
  left_join(
    gene_positions_chr %>% count(chromosome, name = "n_genes"),
    by = "chromosome"
  ) %>%
  left_join(
    gene_clusters %>% count(chromosome, name = "n_clusters"),
    by = "chromosome"
  ) %>%
  left_join(
    superloci %>% count(chromosome, name = "n_superloci"),
    by = "chromosome"
  ) %>%
  mutate(
    n_clusters = replace_na(n_clusters, 0),
    n_superloci = replace_na(n_superloci, 0)
  )

write_tsv(chr_summary, file.path(OUTPUT_DIR, "chromosome_summary.tsv"))

# Family distribution
family_chr_dist <- gene_positions_annotated %>%
  filter(!is.na(family) & !is.na(chromosome)) %>%
  count(chromosome, family, group, name = "n_genes") %>%
  pivot_wider(names_from = family, values_from = n_genes, values_fill = 0)

write_tsv(family_chr_dist, file.path(OUTPUT_DIR, "family_distribution_by_chromosome.tsv"))

# 18) Final summary.
message("\n", paste(rep("=", 70), collapse = ""))
message("ANALYSIS COMPLETE")
message(paste(rep("=", 70), collapse = ""))

message("\nGene Cluster Definition Used:")
message("  - Minimum genes: ", MIN_CLUSTER_GENES)
message("  - Maximum distance: ", MAX_CLUSTER_DISTANCE/1000, " kb")
message("  - Total clusters: ", nrow(gene_clusters))

if (nrow(superloci) > 0) {
  message("\nSuperloci (co-localized multi-family clusters):")
  message("  - Total superloci: ", nrow(superloci))
}

message("\nOutput files in: ", OUTPUT_DIR)
message("\nPlots in: ", PLOTS_DIR)
message("  - 05_chromosome_ideogram_synteny_all_families.pdf (with gene density)")
message("  - 06_chromosome_ideogram_synteny_[group].pdf (with superlocus highlights)")
message("  - 07_gene_clusters_zoom_synteny.pdf (improved labels)")
message("  - 08_chromosome_ideogram_top_DE_genes.pdf")

message("\n", paste(rep("=", 70), collapse = ""))
message("Pipeline completed: ", Sys.time())
message(paste(rep("=", 70), collapse = ""))
