#!/usr/bin/env Rscript
library(argparse)

# parse args
parser <- ArgumentParser("Make single-cell Quality Control Plots")
parser$add_argument(
  "--dir",
  "-d",
  help = "path to run directory",
  default = ""
)
parser$add_argument(
  "--id",
  "-i",
  help = "ID to use for outputs, will read inputs from here"
)
parser$add_argument(
  "--verbose",
  "-v",
  help = "should messages be printed?",
  default = "INFO",
  choices = c("DEBUG", "INFO", "WARN", "ERROR")
)

argv <- parser$parse_args()

# load more libraries
suppressPackageStartupMessages({
  library(here)
  library(future)
  library(Seurat)
  library(patchwork)
  library(tidyverse)
  library(furrr)
  library(log4r)
  library(glue)
})

logger <- logger(threshold = argv$verbose)

source(here("cli/lib/utils.R"))
output_path <- PrepareOutDir(argv)
StopIfOutputDirNotExist(output_path)

plots_path <- here(output_path, "plots")

if (!dir.exists(here(plots_path))) {
  debug(logger, "Plots directory is being created")
  dir.create(here(plots_path))
}

data_filename <- list.files(
  path = here(output_path, "cache/"),
  full.names = TRUE,
  pattern = "^seurat_"
)

if (length(data_filename) > 1) fatal("There is more than one cached object")

debug(logger, paste0("Importing Data from: ", data_filename))
seurat <- readRDS(data_filename)

s <- cc.genes$s.genes
g2m <- cc.genes$g2m.genes

seurat <- CellCycleScoring(seurat, s.features = s, g2m.features = g2m)

cell_cycle1 <- UMAPPlot(seurat, group.by = "Phase")
cell_cycle2 <- UMAPPlot(
  seurat,
  split.by = "Phase",
  combine = FALSE,
  cols = rep_len("black", length(unique(Idents(seurat))))
) + NoLegend()

plot1 <- FeatureScatter(
  seurat,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)
plot2 <- FeatureScatter(
  seurat,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)
std_qc <- plot1 + plot2

hvf_plot <- VariableFeaturePlot(seurat)

pdf(here(glue("{plots_path}/QC_plots.pdf")))
std_qc
hvf_plot
cell_cycle1
cell_cycle2
graphics.off()
