#!/usr/bin/env Rscript
library(argparse)

# parse args
parser <- ArgumentParser("Do GSVA")
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
  "--replicates",
  "-n",
  help = "number of psudoreplicates to use during GSVA",
  default = 3
)
parser$add_argument(
  "--cores",
  "-c",
  help = "number of cores to use",
  default = 5
)
parser$add_argument(
  "--verbose",
  "-v",
  help = "should messages be printed? One of: DEBUG, INFO, WARN, ERROR",
  default = "INFO"
)

argv <- parser$parse_args()

# load more libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(here)
  library(msigdbr)
  library(furrr)
  library(EnhancedVolcano)
  library(tidyverse)
  library(log4r)
})

source("cli/lib/utils.R")
source("cli/lib/gsa.R")

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
diagnosis <- readRDS(data_filename)

# set plan
if (argv$cores == 1) {
  plan("sequential")
} else {
  plan("multisession", workers = argv$cores)
}

gene_sets <- list(
  HALLMARK = msigdbr(species = "Homo sapiens", category = "H"),
  GO = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP"),
  REACTOME = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME"),
  KEGG = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG"),
  ONCO = msigdbr(species = "Homo sapiens", category = "C6")
) %>%
  map(~ ListGeneSets(.x, gene_set = "name", gene_name = "symbol"))

options("future.globals.maxSize" = 10 * 1024^3)

info(logger, "Running GSVA")

furr_options <- furrr_options(
  stdout = FALSE,
  seed = 1000000,
  conditions = structure("condition", exclude = c("message", "warning"))
)

debug(logger, paste0("The Default Assay is: ", DefaultAssay(diagnosis)))

RunGSVAQuietly <- function(...) {
  f <- purrr::quietly(RunGSVA)
  res <- f(...)$result
  return(res)
}
RunStatsQuietly <- function(...) {
  f <- purrr::quietly(RunStats)
  res <- f(...)$result
  return(res)
}
gsva_res <- future_map(
  gene_sets,
  ~ RunGSVAQuietly(
    diagnosis,
    gene_sets = .x,
    replicates = argv$replicates
  ),
  .options = furr_options
)

info(logger, "GSVA Done")

stat_res <- future_map2(
  gsva_res,
  names(gsva_res),
  RunStatsQuietly,
  p = 0.05,
  .options = furr_options
)

file_names <- paste0(names(gene_sets), "_GSVA_volcano_plots.pdf")

top_tables <- map(stat_res, ~ .x[["Top Table"]])

top_tables %>%
  reduce(bind_rows) %>%
  list("All Gene Sets" = .) %>%
  append(top_tables) %>%
  openxlsx::write.xlsx(file = here(output_path, "GSVA_DE_results.xlsx"))

plots <- map(stat_res, ~ .x[["plots"]])

dir.create(here(plots_path, "GSVA"))
walk2(file_names, plots, .f = function(file, plot) {
  pdf(file = here(plots_path, "GSVA", file), onefile = TRUE)
  print(plot)
  graphics.off()
})

source(here("lib/WriteInvocation.R"))
WriteInvocation(argv, output_path = here(output_path, "invocation"))
