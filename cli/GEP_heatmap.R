#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser("Select Genes from CIBERSORTx GEPs for Survival Analysis")
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
  "--verbose",
  "-v",
  help = "should messages be printed? One of: DEBUG, INFO, WARN, ERROR",
  default = "INFO"
)

argv <- parser$parse_args()

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(here)
  library(glue)
  library(log4r)
})

logger <- logger(threshold = argv$verbose)

source(here("cli/lib/utils.R"))
output_path <- PrepareOutDir(argv)
StopIfOutputDirNotExist(output_path)

plots_path <- here(output_path, "plots")

datasets <- c("target_data.txt", "tcga_data.txt", "beat_aml.txt")

ReturnGEP <- function(dataset) {
  fname <- here(
    output_path,
    cibersort_results_dir,
    glue("CIBERSORTxGEP_{dataset}_SM_GEPs_Filtered.txt")
  )
  gep <- readr::read_tsv(fname, col_types = readr::cols())
  return(gep)
}

ht_opt$message <- FALSE

for (set in datasets) {
  do_heatmap <- TRUE
  # Probably not the most efficent
  tryCatch(ReturnGEP(set),
    error = function(e) {
      info(logger, glue("Output for {set} not found"))
      do_heatmap <<- FALSE
    }
  )
  if (do_heatmap) {
    gep <- ReturnGEP(set)
    gep <- column_to_rownames(gep, var = "GeneSymbol")
    gep[is.na(gep)] <- 0
    colnames(gep) <- str_remove(colnames(gep), "cluster")

    col <- circlize::colorRamp2(
      colors = viridis::viridis(n = 3),
      breaks = c(-2, 0, 2)
    )

    gep_sd <- apply(gep, 1, sd)
    q3 <- quantile(gep_sd)[[4]]
    gep_bar <- mean(gep_sd)

    gep_filter <- gep[!(gep_sd <= gep_bar), ]

    gep_scale <- scale(gep_filter)
    pdf(here(plots_path, glue("CIBERSORTxGEP_Heatmap_{set}.pdf")))
    print(Heatmap(gep_scale, col = col, show_row_names = FALSE))
    graphics.off()
  }
}
