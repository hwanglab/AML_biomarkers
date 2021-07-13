#!/usr/bin/env -S Rscript --no-save --quiet
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
  default = "rls_CIBERSORTx"
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

data_filename <- list.files(
  path = here(output_path, "cache"),
  full.names = TRUE,
  pattern = "^deconvoluted_samples_"
)

if (length(data_filename) > 1) fatal("There is more than one cached object")

debug(logger, paste0("Importing Data from: ", data_filename))
deconvoluted_samples <- readRDS(data_filename)

model_use <- argv$model

info(logger, "Preparing Clinical Information")
source("cli/lib/prepare_clin_info.R")

target_deconvoluted <- deconvoluted_samples$TARGET$proportions[[model_use]] %>%
  as_tibble(rownames = "patient_USI") %>%
  separate(patient_USI, into = c(NA, NA, "USI", NA, NA), sep = "-") %>%
  left_join(clinical) %>%
  mutate(across(where(is_character), str_to_lower))

tcga_deconvoluted <- deconvoluted_samples$TCGA$proportions[[model_use]] %>%
  as_tibble(rownames = "case_submitter_id") %>%
  left_join(tcga_ann2) %>%
  filter(flt3_status == "Positive")

beat_aml_decon <- deconvoluted_samples$BeatAML$proportions[[model_use]] %>%
  as_tibble(rownames = "SAMPLE_ID") %>%
  left_join(beat_aml_clinical2) %>%
  filter(FLT3_ITD_CONSENSUS_CALL == "Positive")

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