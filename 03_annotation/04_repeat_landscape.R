#!/usr/bin/env Rscript
#
# Fig. S1b: repeat landscape of O. impluviata - the genome fraction occupied by
# each repeat class, binned by Kimura divergence from its family consensus. Low
# divergence means recent or ongoing activity.
#
# Two versions are drawn, with and without the unclassified interspersed
# repeats, because unclassified elements make up roughly a third of the genome
# and dominate the plot when included.
#
# Run from the repository root:
#   Rscript 03_annotation/04_repeat_landscape.R

# ------------------------------- settings ------------------------------------
DIVSUM_TABLE <- "results/03_annotation/repeats/Oimp.divsum.table"
GENOME_SIZE  <- 474744221      # bp, excluding N/X runs
OUTDIR       <- "results/03_annotation/repeats"
# -----------------------------------------------------------------------------

library(reshape)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(tidyverse)
library(gridExtra)
library(ggpubr)

source("03_annotation/helper.functions.R")   # Terrence Sylvester

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# 1) With unclassified repeats included.
pdf(file.path(OUTDIR, "repeat_landscape_with_unclassified.pdf"), width = 9, height = 5)
plotKimuraDistance(DivsumTab   = DIVSUM_TABLE,
                   genomesSize = GENOME_SIZE,
                   plotUnkowns = TRUE,
                   repClassPlot = TRUE,
                   zoom        = FALSE,
                   savePlot    = FALSE)
dev.off()

# 2) Classified repeats only, which is where the low-divergence peak shows.
pdf(file.path(OUTDIR, "repeat_landscape_no_unclassified.pdf"), width = 9, height = 5)
plotKimuraDistance(DivsumTab   = DIVSUM_TABLE,
                   genomesSize = GENOME_SIZE,
                   plotUnkowns = FALSE,
                   repClassPlot = TRUE,
                   zoom        = FALSE,
                   savePlot    = FALSE)
dev.off()

cat("wrote", file.path(OUTDIR, "repeat_landscape_*.pdf"), "\n")
