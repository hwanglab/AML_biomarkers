#!/usr/bin/env Rscript
source("renv/activate.R")
library(argparser)

# parse args
parser <- arg_parser("Prepare Clinical Data")
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
  help = "ID to use for outputs"
)
parser <- add_argument(
  parser,
  "--cores",
  short = "-c",
  help = "number of cores to use",
  default = 1
)
parser <- add_argument(
  parser,
  "--verbose",
  short = "-v",
  help = "should messages be printed? One of: DEBUG, INFO, WARN, ERROR",
  default = "INFO"
)
parser <- add_argument(
  parser,
  "--model",
  short = "-m",
  help = "which model should be used?",
  default = "CIBERSORTx"
)
parser <- add_argument(
  parser,
  "--batch-corrected",
  short = "-B",
  help = "was batch correction done using CIBERSORTx?",
  flag = TRUE
)

argv <- parse_args(parser)

# load more libraries
suppressPackageStartupMessages({
  library(here)
  library(log4r)
  library(tidyverse)
})

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

logger <- logger(threshold = argv$verbose)

if (argv$model == "CIBERSORTx") {
  info(logger, "Preparing Data from CIBERSORTx")
  data_filename <- c("tcga_data.txt", "target_data.txt", "beat_aml.txt")
  bc_filename <- if_else(argv$batch_corrected, "_Adjusted.txt", "_Results.txt")
  data_filename <- paste0("CIBERSORTx_", data_filename, bc_filename)
  
  data_path <- paste0(output_path, "/cibersort_results/", data_filename)
  debug(logger, paste0("Example Path: ", data_path[[1]]))

  deconvoluted <- data_path %>%
    map(read_tsv, col_types = cols()) %>%
    map(~ select(.x, -`P-value`, -Correlation, -RMSE))
  
  target_deconvoluted <- deconvoluted[[2]]
  tcga_deconvoluted <- deconvoluted[[1]]
  beat_aml_decon <- deconvoluted[[3]]

} else {
  info(logger, "Preparing Data from granulator")
  data_filename <- list.files(
  path = here(output_path, "cache"),
  full.names = TRUE,
  pattern = "^deconvoluted_samples_"
)

if (length(data_filename) > 1) fatal("There is more than one cached object")

debug(logger, paste0("Importing Data from: ", data_filename))
deconvoluted_samples <- readRDS(data_filename)

model_use <- argv$model
target_deconvoluted <- deconvoluted_samples$TARGET$proportions[[model_use]]
tcga_deconvoluted <- deconvoluted_samples$TCGA$proportions[[model_use]]
beat_aml_decon <- deconvoluted_samples$BeatAML$proportions[[model_use]]
}

info(logger, "Preparing Clinical Information")
source("cli/lib/prepare_clin_info.R")

debug(logger, "Done sourcing clinical information tables.")
debug(logger, "Joining TARGET clinical data and deconvolution results")
target_deconvoluted <- target_deconvoluted %>%
  #as_tibble(rownames = "patient_USI") %>%
  separate(Mixture, into = c(NA, NA, "USI", NA, NA), sep = "-") %>%
  left_join(clinical) %>%
  mutate(across(where(is_character), str_to_lower))

debug(logger, "Joining TCGA clinical data and deconvolution results")
tcga_deconvoluted <- tcga_deconvoluted %>%
  #as_tibble(rownames = "case_submitter_id") %>%
  left_join(tcga_ann2, by = c("Mixture" = "Case ID")) %>%
  filter(flt3_status == "Positive")



debug(logger, "Joining BeatAML clinical data and deconvolution results")
beat_aml_decon <- beat_aml_decon %>%
  left_join(beat_aml_clinical2, by = c("Mixture" = "SAMPLE_ID")) %>%
  filter(FLT3_ITD_CONSENSUS_CALL == "Positive")
debug(logger, "Filtering data and saving results")
deconvoluted <- list(
  FLT3 = filter(target_deconvoluted, `FLT3/ITD positive?` == "yes"),
  NEG = filter(
    target_deconvoluted,
    `FLT3/ITD positive?` == "no",
    `CEBPA mutation` == "no"
  ),
  CEBPA = filter(target_deconvoluted, `CEBPA mutation` == "yes"),
  BeatAML = beat_aml_decon,
  TCGA = tcga_deconvoluted
)

saveRDS(deconvoluted, here(output_path, "cache/clinical_deconvoluted.rds"))

source(here("lib/WriteInvocation.R"))
WriteInvocation(argv, output_path = here(output_path, "invocation"))