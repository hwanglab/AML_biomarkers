#!/usr/bin/env Rscript 
source("renv/activate.R")
library(argparse)

# parse args
parser <- ArgumentParser("Prepare reference inputs for CIBERSORTx")
parser$add_argument(
  "--dir",
  "-d",
  help = "path to run directory",
  default = ""
)
parser$add_argument(
  "--id",
  "-i",
  help = "ID to use for outputs"
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
  help = "Should integration be done using Seurat",
  action = "store_true"
)
parser$add_argument(
  "--batch-correct-method",
  "-M",
  "--bc-method",
  help = "which batch correction method to use",
  default = "CCA",
  choices = c("CCA", "RPCA")
)

argv <- parser$parse_args()

# load more libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(here)
  library(xfun)
  library(log4r)
  library(tidyverse)
  library(glue)
})

logger <- logger(argv$verbose)

if (argv$dir == "") {
  output_path <- paste0("outs/", argv$id)
} else {
  output_path <- paste0(parser$run_dir, "/outs/", argv$id)
}
if (!dir.exists(here(output_path))) {
  fatal(logger, "Output directory does not exist")
}

bc_ext <- "no_bc"
if (argv$batch_correct) bc_ext <- argv$batch_correct_method

# read in data
data_filename <- list.files(
  path = here(output_path, "cache/"),
  full.names = TRUE,
  pattern = glue("^{bc_ext}_seurat_dimred_")
)

if (length(data_filename) > 1) fatal("There is more than one cached object")

debug(logger, paste0("Importing Data from: ", data_filename))
diagnosis <- readRDS(data_filename)

seurat_down <- cache_rds(
  {
    subset(diagnosis, downsample = 100)
  },
  file = glue("{bc_ext}_CIBERSORT_ref_prep.rds"),
  dir = paste0(output_path, "/cache/")
)

cibersort_data <- GetAssayData(seurat_down, slot = "data", assay = "SCT")
cibersort_data <- as.data.frame(cibersort_data)
cibersort_clusters <- as.data.frame(
  t(paste0("cluster", Idents(seurat_down)))
)
colnames(cibersort_clusters) <- colnames(cibersort_data)
rownames(cibersort_clusters) <- "GeneSymbol"

cibersort_data2 <- rbind(cibersort_clusters, cibersort_data)

write.table(
  cibersort_data2,
  file = here("cibersort_ref_input.txt"),
  quote = FALSE,
  sep = "\t"
)
system("sed 1d cibersort_ref_input.txt > cibersort_ref_input2.txt")

logical_val <- file.rename(
  here("cibersort_ref_input2.txt"),
  here(output_path, glue("{bc_ext}_cibersort_ref_input.txt"))
)

if (logical_val) debug(logger, "Moving reference succsessful!")

logical_val <- file.remove(here("cibersort_ref_input.txt"))
if (logical_val) debug(logger, "Temporaty Reference File Deleted!")

source(here("lib/WriteInvocation.R"))
WriteInvocation(argv, output_path = here(output_path, "invocation"))
