#!/usr/bin/env Rscript
source("renv/activate.R")
library(argparse)

# parse args
parser <- ArgumentParser("Prepare Clinical Data")
group <- parser$add_mutually_exclusive_group()
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
  "--model",
  "-m",
  help = "which model should be used?",
  default = "CIBERSORTx"
)
parser$add_argument(
  "--batch-correct",
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
group$add_argument(
  "--B-mode",
  "-X",
  help = "was batch correction done using CIBERSORTx?",
  action = "store_true"
)
group$add_argument(
  "--S-mode",
  "-S",
  help = "was batch correction done using CIBERSORTx?",
  action = "store_true"
)

argv <- parser$parse_args()

# load more libraries
suppressPackageStartupMessages({
  library(here)
  library(log4r)
  library(tidyverse)
  library(glue)
  library(scorecard)
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
  bc_filename <- if_else(
    argv$B_mode || argv$S_mode,
    "_Adjusted.txt",
    "_Results.txt"
  )
  data_filename <- paste0("CIBERSORTx_", data_filename, bc_filename)

  if (argv$B_mode) {
    cibersort_results_dir <- "cibersort_results_B"
  } else if (argv$S_mode) {
    cibersort_results_dir <- "cibersort_results_S"
  } else {
    cibersort_results_dir <- "cibersort_results"
  }

  bc_ext <- "no_bc"
  if (argv$batch_correct) bc_ext <- argv$batch_correct_method

  data_path <- glue("{output_path}/{bc_ext}_{cibersort_results_dir}/{data_filename}")
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
  # as_tibble(rownames = "patient_USI") %>%
  separate(Mixture, into = c(NA, NA, "USI", NA, NA), sep = "-") %>%
  left_join(clinical) %>%
  mutate(across(where(is_character), str_to_lower))

debug(logger, "Joining TCGA clinical data and deconvolution results")
tcga_deconvoluted <- tcga_deconvoluted %>%
  # as_tibble(rownames = "case_submitter_id") %>%
  left_join(tcga_ann2, by = c("Mixture" = "Case ID")) %>%
  filter(flt3_status == "Positive")

flt3_target <- filter(target_deconvoluted, `FLT3/ITD positive?` == "yes")
split <- split_df(flt3_target, y = "Event Free Survival Time in Days")

debug(logger, "Joining BeatAML clinical data and deconvolution results")
beat_aml_decon <- beat_aml_decon %>%
  left_join(beat_aml_clinical2, by = c("Mixture" = "SAMPLE_ID")) %>%
  filter(FLT3_ITD_CONSENSUS_CALL == "Positive")
debug(logger, "Filtering data and saving results")
deconvoluted <- list(
  TRAIN = split$train,
  FLT3 = split$test,
  NEG = filter(
    target_deconvoluted,
    `FLT3/ITD positive?` == "no",
    `CEBPA mutation` == "no"
  ),
  CEBPA = filter(target_deconvoluted, `CEBPA mutation` == "yes"),
  BeatAML = beat_aml_decon,
  TCGA = tcga_deconvoluted
)

saveRDS(deconvoluted, here(output_path, glue("cache/{bc_ext}_clinical_deconvoluted.rds")))

source(here("lib/WriteInvocation.R"))
WriteInvocation(argv, output_path = here(output_path, "invocation"))
