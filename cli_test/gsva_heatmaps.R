#!/usr/bin/env Rscript
source("renv/activate.R")
library(argparse)

# parse args
parser <- ArgumentParser("Make GSVA Heatmaps")
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
  "--clusters",
  "-c",
  help = "what cluster(s) to plot",
  nargs = "+"
)
parser$add_argument(
  "--filter-genesets",
  "-f",
  help = "should the gene sets be filtered",
  action = "store_true"
)
parser$add_argument(
  "--p-value",
  "-p",
  help = "p value to use for filtering",
  default = 0.01
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
  library(here)
  suppressWarnings(library(rsinglecell))
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

if (!any(is.na(argv$clusters))) {
  info(logger, "Selecting Clusters")
  data_filtered <- map(data, filter, cluster %in% argv$clusters)
  file_name_ext <- paste0(argv$clusters, collapse = "_")
} else {
  data_filtered <- data
  file_name_ext <- "all_clusters"
}

GetSignificantGeneSets <- function(df, p = 0.01) {
  data <- filter(df, adj.P.Val <= p)
  res <- unique(data$geneset)
  return(res)
}

if (argv$filter_genesets) {
  info(logger, "Selecting Gene Sets")

  # p <- if_else(is.na(argv$p_value), 0.01, argv$p_value, missing = 0.01)
  p <- argv$p_value
  debug(logger, paste0("--p-value is parsed as: ", argv$p_value, "."))
  info(logger, paste0("p is set to ", p, "."))
  pos_is_na <- possibly(is.na, otherwise = TRUE)
  if (pos_is_na(p) || length(p) == 0) {
    error(logger, "The p value cutoff is somehow not set")
    quit(status = 1)
  }
  sets <- map(data_filtered, GetSignificantGeneSets, p = p)
  nsets <- map_int(sets, length)
  debug(logger, paste0("The number of sets is between ", min(nsets), " and ", max(nsets)))
  data_filtered <- map2(data_filtered, sets, ~ filter(.x, geneset %in% .y))
  file_name_ext <- paste0(file_name_ext, "filtered_p=", p)
}

data_filtered <- data_filtered %>%
  map(select, geneset, logFC, cluster) %>%
  map(pivot_wider, names_from = cluster, values_from = logFC) %>%
  map(column_to_rownames, var = "geneset") %>%
  map(as.matrix)

plots <- list(0)

cluster_cols <- if_else(length(argv$clusters) == 1, TRUE, FALSE)

for (i in seq_along(data_filtered)) {
  filename <- paste0(names(data_filtered)[[i]], "_heatmap_(", file_name_ext, ").pdf")
  plots[[i]] <- pheatmap(
    data_filtered[[i]],
    color = viridis::magma(n = 100),
    fontsize_row = 5,
    cellwidth = 10,
    filename = here(plots_path, filename),
    cluster_rows = !argv$filter_genesets,
    cluster_cols = cluster_cols
  )
}

source(here("lib/WriteInvocation.R"))
WriteInvocation(argv, output_path = here(output_path, "invocation"))
