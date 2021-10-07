#!/usr/bin/env Rscript
source("renv/activate.R")
library(argparse)

# parse args
parser <- ArgumentParser("Assess Single Cell Gene Correlation")
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
  help = "should messages be printed? One of: DEBUG, INFO, WARN, ERROR",
  default = "INFO"
)
parser$add_argument(
  "--batch-correct",
  "-X",
  help = "Should intetegraed data be used from Seurat",
  action = "store_true"
)
parser$add_argument(
  "--genes",
  "-g",
  help = "genes to assess correlation of",
  nargs = "+"
)
parser$add_argument(
  "--gene-et-name",
  "-G",
  help = "Name of gene set for filenames"
)
suppressPackageStartupMessages({
  library(Seurat)
  library(here)
  library(pheatmap)
  library(viridis)
  library(log4r)
  library(tidyverse)
  library(glue)
})

logger <- logger(threshold = argv$verbose)

if (argv$batch_correct) {
data_filename <- list.files(
  path = here(output_path, "cache/"),
  full.names = TRUE,
  pattern = "^seurat_integrated_bc_"
)
} else {
  data_filename <- list.files(
  path = here(output_path, "cache/"),
  full.names = TRUE,
  pattern = "^seurat_dimred_"
)
}
if (length(data_filename) > 1) {
  fatal(logger, "There is more than one cached object")
  quit()
}

create.dir(glue("{plots_path}/leading_edge_{argv$gene_set_name}"))
seurat <- readRDS(data_filename)

info(logger, "Averaging Expression")

expression <- AverageExpression(seurat, features = argv$genes, assays = "SCT", slot = "data", group.by = "prognosis", verbose = FALSE)$SCT
filename_feats <- glue_collapse(features, sep = "_", width = 10)
filename <- glue("{plots_path}/leading_edge_{argv$gene_set_name}/heatmap_{filename_feats}.tsv")

info(logger, "Plotting Results")
pheatmap(expression, color = viridis(100), cellwidth = 10, cellheight = 10, cluster_rows = FALSE, cluster_cols = FALSE, filename = filename)