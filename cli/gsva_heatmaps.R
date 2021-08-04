#!/usr/bin/env Rscript
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
  "--clusters",
  short = "-c",
  help = "what cluster(s) to plot",
  nargs = Inf
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
  library(here)
  library(rsinglecell)
  library(pheatmap)
  library(readxl)
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
if (!dir.exists(here(output_path))) {
  fatal(logger, "Output directory does not exist")
}

# read in data
data_filename <- here(output_path, "GSVA_DE_results.xlsx")
debug(logger, paste0("Importing Data from: ", data_filename))
sheet_names <- excel_sheets(data_filename) %>% str_subset("[^All Gene Sets]")

data <- map(sheet_names, ~ as.data.frame(read_excel(data_filename, sheet = .x)))
names(data) <- sheet_names

data_filtered <- map(data, filter, cluster %in% argv$clusters)

data_filtered <- data_filtered %>%
  map(select, geneset, logFC, cluster) %>%
  map(pivot_wider, names_from = cluster, values_from = logFC) %>%
  map(column_to_rownames, var = "geneset") %>%
  map(as.matrix)

plots <- list(0)

for (i in seq_along(data_filtered)) {
    filename <- paste0(names(data_filtered)[[i]], "_heatmap.pdf")
  plots[[i]] <- pheatmap(
    data_filtered[[i]],
    color = viridis::magma(n = 100),
    fontsize_row = 5,
    cellwidth = 10,
    filename = here(plots_path, filename)
  )
}

source(here("lib/WriteInvocation.R"))
WriteInvocation(argv, output_path = here(output_path, "invocation"))
