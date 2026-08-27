#!/usr/bin/env Rscript
#
# Orthology distribution figure: the nine orthogroup categories from step 08,
# stacked per species and ordered to match the tips of the species tree.
#
# Run from the repository root:
#   Rscript 04_orthology/10_orthology_histogram.R

# ------------------------------- settings ------------------------------------
PROJECT     <- Sys.getenv("PROJECT", "results/04_orthology")
PLOTS_DIR   <- file.path(PROJECT, "plots")
CSV_WIDE    <- file.path(PLOTS_DIR, "orthology_histogram_counts_wide.csv")
ORDER_TXT   <- file.path(PLOTS_DIR, "orthology_tip_order.txt")
                       # tip order from step 08, NOT species_taxonomy.tsv
OUTPUT_PLOT <- "orthology_histogram.pdf"
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

PLOT_WIDTH <- 8
PLOT_HEIGHT <- 5

# 1) The nine orthology categories, in plotting order, with their colours.
cat_levels <- c(

  "SingleCopy_Universal",

  "Multicopy_Universal",

  "Species_Specific",

  "Curculionoidea_Specific",

  "Cerambycidae_Specific",

  "Chrysomelidae_Specific",

  "Phytophaga_Specific",

  "Other_Orthologs",

  "Unassigned_Genes"

)

cat_labels <- c(

  "1:1:1",

  "N:N:N",

  "Species-specific",

  "Curculionoidea-specific",

  "Cerambycidae-specific",

  "Chrysomelidae-specific",

  "Phytophaga-specific",

  "Other orthologs",

  "Unassigned genes"

)

cat_colors <- c(

  "SingleCopy_Universal" =
    "#1f78b4",

  "Multicopy_Universal" =
    "#e31a1c",

  "Species_Specific" =
    "#33a02c",

  "Curculionoidea_Specific" =
    "#b15928",

  "Cerambycidae_Specific" =
    "#6a3d9a",

  "Chrysomelidae_Specific" =
    "#ff7f00",

  "Phytophaga_Specific" =
    "#a6cee3",

  "Other_Orthologs" =
    "#C7C7C7",

  "Unassigned_Genes" =
    "#666666"

)

# Helpers for labelling and ordering.
clean_species_name <- function(x) {

  x <- gsub(
    "_protein.*",
    "",
    x
  )

  x <- gsub(
    "_",
    " ",
    x
  )

  return(x)

}

# 2) Load the wide count table and the tip order from step 08.
message(
  paste(
    rep("=",70),
    collapse=""
  )
)

message(
  "Orthology Distribution Analysis"
)

message(
  paste(
    rep("=",70),
    collapse=""
  )
)

setwd(WORK_DIR)

message("Loading input files...")

message(
  "  Data: ",
  CSV_WIDE
)

message(
  "  Order: ",
  ORDER_TXT
)

df_wide <- read_csv(
  CSV_WIDE,
  show_col_types = FALSE
)

tip_order <- read_lines(
  ORDER_TXT
)

tip_order <- tip_order[
  tip_order != ""
]

message(
  "Species in CSV: ",
  nrow(df_wide)
)

message(
  "Species in order file: ",
  length(tip_order)
)

# 3) Check that the table and the tip order name the same species.
missing_species <- setdiff(
  df_wide$Species,
  tip_order
)

if(length(missing_species) > 0){

  warning(
    paste(
      "Species in CSV not found in order file:",
      paste(
        missing_species,
        collapse=", "
      )
    )
  )

}

# 4) Reshape to long form, one row per species and category.
df_long <- df_wide %>%

  pivot_longer(

    cols = -Species,

    names_to = "Category",

    values_to = "Count"

  )

# Match exact species order

df_long$Species <- factor(

  df_long$Species,

  levels = rev(tip_order)

)

df_long$Pretty <- factor(

  clean_species_name(
    as.character(df_long$Species)
  ),

  levels = rev(
    clean_species_name(tip_order)
  )

)

df_long$Category <- factor(

  df_long$Category,

  levels = cat_levels

)

message(
  "Data prepared for plotting."
)

# 5) Draw the stacked histogram.
message(
  "Generating plot..."
)

p <- ggplot(

  df_long,

  aes(
    x = Count,
    y = Pretty,
    fill = Category
  )

) +

geom_col(

  width = 0.7,

  position = position_stack(
    reverse = TRUE
  )

) +

scale_fill_manual(

  values = cat_colors,

  breaks = cat_levels,

  labels = cat_labels,

  name = "Orthology category"

) +

theme_classic() +

theme(

  axis.text.y =
    element_text(
      face="italic",
      size=9
    ),

  axis.text.x =
    element_text(
      size=9
    ),

  legend.text =
    element_text(
      size=8
    ),

  legend.title =
    element_text(
      size=9
    ),

  legend.position="right"

) +

labs(

  x="Number of genes",

  y=NULL

)

output_path <- file.path(
  PLOTS_DIR,
  OUTPUT_PLOT
)

ggsave(

  output_path,

  p,

  width=PLOT_WIDTH,

  height=PLOT_HEIGHT,

  dpi=300

)

message(
  "Plot saved:"
)

message(
  output_path
)

# 6) Report what was written.
summary <- df_long %>%

  group_by(Category) %>%

  summarise(

    total_genes=sum(Count),

    mean_genes=round(mean(Count),1),

    median_genes=median(Count),

    .groups="drop"

  )

print(summary)

message(
  "Analysis complete."
)
