#!/usr/bin/env Rscript
library(argparse)

# parse args
parser <- ArgumentParser("Prepare Clinical Data")
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
  "--additional-deconvolutions",
  "-A",
  metavar = "BASENAME",
  help = glue::glue("Additional clinical_deconvolution.rds files from \\
  external runs of this script. Useful for PHI from collaborators. These \\
  files should be placed in outs/ID/external. No need for an rds extension. \\
  You can also provide a full path")
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
} else {
  output_path <- paste0(parser$run_dir, "/outs/", argv$id)
}
if (!dir.exists(here(output_path))) {
  fatal(logger, "Output directory does not exist")
}

logger <- logger(threshold = argv$verbose)

if (argv$model == "CIBERSORTx") {
  info(logger, "Preparing Data from CIBERSORTx")
  data_filename <- c("tcga_data.txt", "target_data.txt", "beat_aml.txt")

  is_batch_corrected <- file.exists(here(glue("{output_path}/cibersort_results/b_mode"))) || file.exists(here(glue("{output_path}/s_mode")))

  debug(logger, glue("Detected data that is {if_else(is_batch_corrected, '', 'not ')}batch corrected"))

  bc_filename <- if_else(
    is_batch_corrected,
    "_Adjusted.txt",
    "_Results.txt"
  )
  data_filename <- paste0("CIBERSORTx_", data_filename, bc_filename)

  data_path <- glue("{output_path}/cibersort_results/{data_filename}")
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
  left_join(clinical, by = "USI") %>%
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
  "TARGET:FLT3" = split$test,
  "TARGET:NEG" = filter(
    target_deconvoluted,
    `FLT3/ITD positive?` == "no",
    `CEBPA mutation` == "no"
  ),
  "TARGET:CEBPA" = filter(target_deconvoluted, `CEBPA mutation` == "yes"),
  "BeatAML:FLT3" = beat_aml_decon,
  "TCGA:FLT3" = tcga_deconvoluted
)

adtnl_dec <- argv$additional_deconvolutions
external_data <- list()

TestIfPath <- function(path) {
  tryCatch(
    {
      res <- readRDS(path)
      info(logger, "Successfully read file from path")
      return(res)
    },
    error = function(e) {
      warn(logger, "File not valid system path.")
      info(
        logger,
        glue("Did you put the file in {here(output_path, 'external')}/ or provide the full path?")
      )
    }
  )
}

if (!is.null(adtnl_dec)) {
  debug(
    logger,
    glue("Adding in additional deconvolutions: {glue_collapse(adtnl_dec)}")
  )
  for (i in seq_along(adtnl_dec)) {
    info(logger, glue("Loading {adtnl_dec[[i]]}"))
    suppressWarnings({
      tryCatch(
        external_data[[i]] <- readRDS(
          here(output_path, "external", glue("{adtnl_dec[[i]]}.rds"))
        ),
        error = function(e) {
          info(
            logger,
            glue("File not found. Testing if full path.")
          )
          external_data[[i]] <- TestIfPath(adtnl_dec[[i]])
        }
      )
    })
  }
  if (!is_empty(external_data)) {
    ext_dat <- unlist(external_data, recursive = FALSE)
    deconvoluted <- append(deconvoluted, ext_dat)
  }
}

saveRDS(deconvoluted, here(output_path, glue("cache/clinical_deconvoluted.rds")))

source(here("lib/WriteInvocation.R"))
WriteInvocation(argv, output_path = here(output_path, "invocation"))
