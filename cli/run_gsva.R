#!/usr/bin/env -S Rscript --no-save --quiet
source("renv/activate.R")
library(argparser)

# parse args
parser <- arg_parser("Do GSVA")
parser <- add_argument(
  parser, "--dir",
  short = "-d",
  help = "path to run directory",
  default = ""
)
parser <- add_argument(
  parser,
  "--id",
  short = "-i",
  help = "ID to use for outputs, will read inputs from here"
)
parser <- add_argument(
  parser,
  "--replicates",
  "-n",
  help = "number of psudoreplicates to use during GSVA",
  default = 3
)
parser <- add_argument(
  parser,
  "--cores",
  "-c",
  help = "number of cores to use",
  default = 5
)
parser <- add_argument(
  parser,
  "--verbose",
  "-v",
  help = "should messages be printed? One of: DEBUG, INFO, WARN, ERROR",
  default = "INFO"
)

argv <- parse_args(parser)

# load more libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(here)
  suppressWarnings(
    library(rsinglecell)
  )
  library(msigdbr)
  library(furrr)
  library(EnhancedVolcano)
  library(tidyverse)
  library(log4r)
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
  seed = 1000000
)

gsva_res <- future_map(
  gene_sets,
  ~ RunGSVA(
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
  RunStats,
  p = 0.05,
  .options = furr_options
)

file_names <- paste0(names(gene_sets), "_GSVA_volcano_plots.pdf")

for (i in seq_along(stat_res)) {
  pdf(file = here(plots_path, file_names[[i]]), onefile = TRUE)
  stat_res[[i]]$plots
  graphics.off()
}

top_tables <- map(stat_res, ~ .x[["Top Table"]])

top_tables %>%
  reduce(bind_rows) %>%
  list("All Gene Sets" = .) %>%
  append(top_tables) %>%
  openxlsx::write.xlsx(file = here(output_path, "GSVA_DE_results.xlsx"))

source(here("lib/WriteInvocation.R"))
WriteInvocation(argv, output_path = here(output_path, "invocation"))
