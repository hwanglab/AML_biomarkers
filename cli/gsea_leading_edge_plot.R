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
  "--gene-set-name",
  "-G",
  help = "Name of gene set for filenames"
)
parser$add_argument(
  "--cluster",
  "-c",
  help = "number of cluster to use"
)
parser$add_argument(
  "--dataset",
  "-D",
  help = "which dataset to use",
  required = TRUE
)
argv <- parser$parse_args()
suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(Seurat)
  library(here)
  library(pheatmap)
  library(viridis)
  library(log4r)
  library(tidyverse)
  library(glue)
})

logger <- logger(threshold = argv$verbose)

if (argv$dir == "") {
  output_path <- paste0("outs/", argv$id)
  plots_path <- paste0("plots/", argv$id)
} else {
  output_path <- paste0(parser$run_dir, "/outs/", argv$id)
  plots_path <- paste0(parser$run_dir, "/plots/", argv$id)
}

bc_ext <- ""
if (argv$batch_correct) bc_ext <- "bc_"

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

### new method
# Get GSEA Res
# Get Leading Edge of Pathway from contrast
# Translate from entrez ID (if applicable)
# Use genes as input to Heatmap

debug(logger, "Reading in GSVA results")
gsea_res <- read_tsv(glue("{output_path}/{bc_ext}GSEA.tsv"), col_types = cols())

gsea_res <- gsea_res %>% filter(ID == argv$gene_set_name)

if (nrow(gsea_res) == 1) {
  fatal(logger, glue("There are no rows matching {argv$gene_set_name}"))
  quit(staus = 1)
}

debug(logger, "Extracting Leading Edge Genes")
gsea_res <- gsea_res %>% filter(prognosis == argv$dataset)

# if (length(argv$cluster) != 0) {
#   gsea_res <- filter(gsea_res, cluster == argv$cluster)
# }
if (nrow(gsea_res) == 0) {
  fatal(
    logger,
    glue("There are no rows matching {argv$dataset} and/or {argv$cluster}")
  )
  quit(status = 1)
}
leading_edge <- pull(gsea_res, core_enrichment)

leading_edge <- unlist(strsplit(leading_edge, "/"))

database <- pull(gsea_res, gene_set)

if (database != "GO") {
  debug(logger, "Translating entrez IDs to symbols")
  leading_edge <- mapIds(org.Hs.eg.db, leading_edge, "SYMBOL", "ENTREZID")
}

debug(logger, glue("Saving to: {plots_path}/leading_edge_{argv$gene_set_name}"))
dir.create(
  glue("{plots_path}/leading_edge/{argv$dataset}"),
  recursive = TRUE,
  showWarnings = FALSE
)
seurat <- readRDS(data_filename)

info(logger, "Averaging Expression")

expression <- AverageExpression(
  seurat,
  features = leading_edge,
  assays = "SCT",
  slot = "data",
  group.by = "prognosis",
  verbose = FALSE
)$SCT
filename_feats <- glue_collapse(leading_edge, sep = "_", width = 10)
filename <- glue(
  "{plots_path}/leading_edge/{argv$dataset}/heatmap_{argv$gene_set_name}.pdf"
)

cellheight <- if_else(length(leading_edge) > 10, 1, 10)
info(logger, "Plotting Results")
pheatmap(
  expression,
  color = viridis(100),
  cellwidth = 10,
  cellheight = cellheight,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_row = 1,
  filename = filename
)
info(logger, "Done!")