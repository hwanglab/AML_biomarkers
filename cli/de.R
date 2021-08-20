#!/usr/bin/env Rscript
source("renv/activate.R")
library(argparse)

# parse args
parser <- ArgumentParser("Do Differential Expression")
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
parser$add_argument(
  "--cores",
  "-c",
  help = "how many cores should be used",
  default = 1
)
parser$add_argument(
  "--clusters",
  "-C",
  help = "clusters to use for DE",
  nargs = "+"
)
parser$add_argument(
  "--sig-level",
  "-p",
  help = "p-value threshold to use",
  default = 0.05,
  metavar = "P"
)

argv <- parser$parse_args()

# load more libraries
suppressPackageStartupMessages({
  library(here)
  library(future)
  library(Seurat)
  suppressWarnings(library(rsinglecell))
  library(pheatmap)
  library(readxl)
  library(tidyverse)
  library(furrr)
  library(log4r)
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

if (!dir.exists(here(plots_path))) {
  debug(logger, "Plots directory is being created")
  dir.create(here(plots_path))
}

if (!dir.exists(here(output_path))) {
  fatal(logger, "Output directory does not exist")
}

# read in data
data_filename <- list.files(
  path = here(output_path, "cache/"),
  full.names = TRUE,
  pattern = "^seurat_"
)

if (length(data_filename) > 1) fatal("There is more than one cached object")

debug(logger, paste0("Importing Data from: ", data_filename))
seurat <- readRDS(data_filename)

# set plan
if (argv$cores == 1) {
  plan("sequential")
} else {
  plan(tweak("multicore", workers = as.integer(argv$cores)))
}

info(logger, "Printing Plots")
WilcoxTestPossibly <- purrr::possibly(wilcox.test, otherwise = list(p.value = NA))

FindClusterFreq <- function(object, metadata, cluster, sort_by = NULL, .debug = FALSE) {
  # quote expressions
  md <- rlang::enquos(metadata)
  c <- rlang::enquo(cluster)
  if (.debug == TRUE) browser()
  # select meta.data columns from seurat object and calculate freq
  if (is.null(sort_by)) {
    freq <- object %>%
      dplyr::group_by(dplyr::across(.cols = all_of(c(metadata, cluster)))) %>%
      summarise(n = n(), .groups = "keep") %>%
      mutate(freq = n / sum(n) * 100)
  } else {
    sb <- rlang::enquo(sort_by)
    freq <- object %>%
      dplyr::group_by(dplyr::across(.cols = all_of(c(metadata, cluster)))) %>%
      summarise(n = n(), .groups = "keep") %>%
      mutate(freq = n / sum(n) * 100) %>%
      dplyr::arrange(!!sb)
  }
  return(freq)
}
freq <- FindClusterFreq(
  seurat[[]],
  c("patient_id", "prognosis"),
  "clusters"
)

freq_stat <- freq %>%
  group_by(clusters) %>%
  summarise(
    pval = WilcoxTestPossibly(
      x = freq,
      y = as.numeric(as.factor(prognosis))
    )$p.value, .groups = "keep"
  )

wilcox_clusters <- freq_stat %>%
  filter(pval <= argv$sig_level) %>%
  pull(clusters)

print_clusters <- glue_collapse(wilcox_clusters, sep = ", ", last = " and ")
debug(logger, glue("The following clusters are selected: {print_clusters}"))

bad_cells <- WhichCells(seurat, idents = wilcox_clusters, invert = TRUE)
Idents(seurat, cells = bad_cells) <- "Not Significant"

pdf(here(plots_path, "UMAP.pdf"), onefile = TRUE, width = 10)
DimPlot(seurat, reduction = "umap", group.by = "clusters")
DimPlot(seurat, reduction = "umap", split.by = "prognosis")
graphics.off()

pdf(here(plots_path, "feature_heatmap.pdf"), width = 18)
DoHeatmap(seurat, group.by = "clusters")
graphics.off()

if (!file.exists(here(output_path, "cluster_differential_expression.tsv"))) {
  info(loggger, "Doing differential Expression")
  markers <- FindAllMarkers(seurat, test.use = "MAST")
  write_tsv(markers, file = here(output_path, "cluster_differential_expression.tsv"))
}

if (!is.na(argv$clusters)) {
  if (length(argv$clusters) == 1) {
    error(logger, "You must have at least 2 clusters selected")
    quit()
  }
  cells <- WhichCells(seurat, idents = argv$clusters)
  seurat_selected <- seurat[, cells]
  markers <- FindAllMarkers(seurat_selected, test.use = "MAST")
  clusters <- glue_collapse(clusters, sep = ", ", final = " and ")
  fname <- glue("cluster_differential_expression_({clusters}).tsv")
  write_tsv(markers, file = here(output_path, fname))
}
