#!/usr/bin/env -S Rscript --no-save --quiet
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
  "--related-test-id",
  "-R",
  help = "ID used in previous run to pull lasso model information from",
  metavar = "ID"
)
parser$add_argument(
  "--test-id",
  "-I",
  help = "ID to use for tests",
  default = "incremental",
  metavar = "ID"
)
parser$add_argument(
  "--verbose",
  "-v",
  help = "should messages be printed?",
  default = "INFO",
  choices = c("DEBUG", "INFO", "WARN", "ERROR")
)
parser$add_argument(
  "--training-data",
  "-t",
  help = "what should be used to train on?",
  default = "FLT3",
  choices = c("FLT3", "CEBPA", "NEG", "BeatAML", "TCGA")
)
parser$add_argument(
  "--training-split",
  "-T",
  help = "should the training data be split? 0 = no split",
  default = 0.7,
  metavar = "split"
)
parser$add_argument(
  "--training-axis",
  "-a",
  help = "What time axis to use for training",
  default = "Event Free Survival Time in Days",
  metavar = "axis"
)
parser$add_argument(
  "--training-vars",
  "-V",
  help = "What vars to include during Training. Should be a quoted R expression",
  metavar = "vars"
)
parser$add_argument(
  "--exclude",
  "-E",
  help = "Should training-vars be excluded instead?",
  action = "store_true"
)
parser$add_argument(
  "--split",
  "-S",
  help = "Should previous split information be used",
  action = "store_true"
)
parser$add_argument(
  "--datasets",
  "-D",
  help = "Which datasets should be included in addtion to TARGET?",
  choices = c("BeatAML", "TCGA"),
  nargs = "+"
)

argv <- parser$parse_args()

suppressPackageStartupMessages({
  library(here)
  library(log4r)
  library(survival)
  library(survminer)
  library(scorecard)
  library(glmnet)
  suppressWarnings(library(rsinglecell))
  library(tidyverse)
  library(glue)
})

logger <- logger(threshold = argv$verbose)

if (argv$dir == "") {
  output_path <- paste0("outs/", argv$id)
  output_path_root <- output_path
} else {
  output_path <- paste0(parser$run_dir, "/outs/", argv$id)
  output_path_root <- output_path
}

plots_path <- here(output_path, "plots")

if (argv$test_id == "incremental") {
  dir_names <- list.dirs(here(output_path), full.names = FALSE)
  debug(
    logger,
    paste0(
      "The following directories exist: ",
      paste0(dir_names, collapse = ", ")
    )
  )
  dir_names <- as.numeric(dir_names)
  dir_names <- dir_names[!is.na(dir_names)]
  if (is_empty(dir_names)) {
    dir_names <- 0
  }

  debug(
    logger,
    paste0(
      "The following directories exist: ",
      paste0(dir_names, collapse = ", ")
    )
  )
  new_name <- max(dir_names, na.rm = TRUE) + 1
  info(logger, c("Using ", new_name, " as the test id"))
  output_path <- here(output_path, new_name)
  plots_path <- here(plots_path, new_name)
  dir.create(output_path, showWarnings = FALSE)
  dir.create(plots_path, recursive = TRUE, showWarnings = FALSE)
} else {
  output_path <- here(output_path, argv$test_id)
  plots_path <- here(plots_path, argv$test_id)
  info(logger, c("Using ", argv$test_id, " as the test id"))
  dir.create(output_path, showWarnings = FALSE)
  dir.create(plots_path, recursive = TRUE, showWarnings = FALSE)
}

if (!dir.exists(here(output_path))) {
  fatal(logger, "Output directory does not exist")
}

dataset_names <- c("target", "tcga", "beat_aml")

datasets <- map(
  dataset_names,
  ~ read_tsv(
    glue("cibersort_in/{.x}_data.txt"),
    col_types = cols()
  )
)
names(datasets) <- dataset_names

cibersort_gep <- read_tsv(
  glue(
    "{output_path_root}/cibersort_results/CIBERSORTx_cell_type_sourceGEP.txt"
  ),
  col_types = cols()
)

PrintGluedList <- function(l) glue_collapse(l, sep = ", ", last = " and ", width = 20)

if (!is.na(argv$related_test_id)) {
  old_lasso_model <- read_tsv(
    glue("{output_path_root}/{argv$related_test_id}/lasso_model_coefs.tsv"),
    col_types = cols()
  )
  old_lasso_model <- filter(old_lasso_model, s0 != 0)
  clusters <- str_subset(old_lasso_model$rownames, "^cluster")
  clusters_print <- str_remove(clusters, "cluster")
  info(logger, glue("Using the following clusters: {PrintGluedList(clusters_print)}"))
} else {
  clusters <- paste0("cluster", argv$clusters)
}

cibersort_gep_filtered <- select(
  cibersort_gep,
  all_of(clusters),
  genesymbols
) %>%
  mutate(row_sum = rowSums(select(., -genesymbols))) %>%
  filter(row_sum != length(clusters)) %>%
  select(-row_sum)

genes <- pull(cibersort_gep_filtered, genesymbols)

info(logger, "Pulling Genes of Interest from datasets")

# suppress: the `x` argument of `as_tibble.matrix()` must have unique column
# names if `.name_repair` is omitted as of tibble 2.0.0.
# Using compatibility `.name_repair`.
suppressWarnings({
  datasets <- map(datasets, ~ filter(.x, GeneSymbol %in% genes)) %>%
    map(transposeDF)
})

source("cli/lib/prepare_clin_info.R")

datasets[["target"]]$samples <- datasets[["target"]]$samples %>%
  str_remove("TARGET-[:digit:]+-") %>%
  str_remove("-[:alnum:]{3}-[:alnum:]{3}")

clinical <- rename(clinical, samples = USI) %>%
  mutate(
    across(
      .cols = c("FLT3/ITD positive?", "CEBPA mutation"),
      .fns = tolower
    )
  )
tcga_ann <- rename(tcga_ann, samples = case_submitter_id)
beat_aml_clinical <- rename(beat_aml_clinical, samples = SAMPLE_ID)

clinical_list <- list(clinical, tcga_ann, beat_aml_clinical)
annotated_datasets <- map2(datasets, clinical_list, inner_join, by = "samples")

debug(logger, glue("TARGET has {nrow(annotated_datasets[['target']])} rows and {ncol(annotated_datasets[['target']])} cols"))

target_subsets <- list(
  FLT3 = filter(annotated_datasets[["target"]], `FLT3/ITD positive?` == "yes"),
  NEG = filter(
    annotated_datasets[["target"]],
    `FLT3/ITD positive?` == "no",
    `CEBPA mutation` == "no"
  ),
  CEBPA = filter(annotated_datasets[["target"]], `CEBPA mutation` == "yes")
)

if (any(map_lgl(target_subsets, plyr::empty))) {
  error(logger, "TARGET subseting FAILED")
  quit(status = 1)
}

debug(logger, glue("Names of initial data sets: {PrintGluedList(names(annotated_datasets))}"))

names(annotated_datasets) <- c("TARGET", "TCGA", "BeatAML")

data <- target_subsets

if ("TCGA" %in% argv$datasets) {
  data <- append(
    annotated_datasets[2],
    data
  )
}
if ("BeatAML" %in% argv$datasets) {
  data <- append(
    annotated_datasets[3],
    data
  )
}

debug(logger, glue("Names of final data sets: {PrintGluedList(names(data))}"))

# info(logger, "Preparing scaled data using CIBERSORTx GEP")

# goal is to get the row means for the clusters and multiply those with the datasets

# cibersort_gep_filtered[cibersort_gep_filtered == 1] <- NA

# cibersort_gep_filtered <- cibersort_gep_filtered %>%
#   mutate(
#     row_mean = rowSums(
#       select(., -genesymbols),
#       na.rm = TRUE
#     ) / (ncol(cibersort_gep_filtered) - 1)
#   )
# cibersort_gep_filtered <- select(cibersort_gep_filtered, -matches("cluster"))

# temp_data <- data$FLT3
# res <- list(0)

# CIBERSORTxScaleFactor <- function(data, gene, gep) {
#   dat <- as_vector(data[, gene])
#   scale_factor <- gep[which(gep$genesymbols == gene), "row_mean"] %>%
#   pull(row_mean)
#   return(dat * scale_factor)
# }
# PossiblyCIBERSORTxScaleFactor <- possibly(CIBERSORTxScaleFactor, otherwise = NULL)

# map(
#   cibersort_gep_filtered$genesymbols,
#   ~ PossiblyCIBERSORTxScaleFactor(temp_data, .x, cibersort_gep_filtered)
# ) %>% print()

# TODO: Row Bind

debug(logger, "Splitting Data")
if (argv$split &
  !file.exists(
    here("outs", argv$id, "cache", paste0(argv$training_data, "_rownames.rds"))
  )
) {
  warn(logger, "The argument split was used, but a cached object was not found. Ignoring option")
  argv$split <- FALSE
}

if (argv$split) {
  debug("Reading files for split")
  rn_split <- readRDS(here("outs", argv$id, "cache", paste0(argv$training_data, "_rownames.rds")))
  split <- map(rn_split, ~ data[[argv$training_data]][
    data[[argv$training_data]][[1]] %in% .x,
  ])
  names(split) <- c("train", "test")
} else {
  debug(
    logger,
    paste0(
      "Using ratios: ",
      argv$training_split,
      " and ",
      1 - argv$training_split
    )
  )
  split <- split_df(
    data[[argv$training_data]],
    y = argv$training_axis,
    ratios = c(argv$training_split, 1 - argv$training_split)
  )
  if (plyr::empty(split$train)) {
    error(logger, "Split Failed. Training Data empty")
    quit(status = 1)
  }
  saveRDS(
    map(split, ~ .x[[1]]),
    file = here("outs", argv$id, "cache", paste0(argv$training_data, "_rownames.rds"))
  )
}

debug(logger, "Preparing Data for training")
debug(logger, glue("Using {argv$training_data} as the training data"))

if (argv$training_split != 0) {
  training <- split$train
  data[[argv$training_data]] <- NULL
  data[[paste0(argv$training_data, ": Training")]] <- split$train
  data[[paste0(argv$training_data, ": Test")]] <- split$test
  debug(logger, glue("Training Data: {nrow(split$training)} rows and {ncol(split$training)} cols"))
  debug(logger, glue("Test Data: {nrow(split$test)} rows and {ncol(split$test)} cols"))
} else {
  training <- data[[argv$training_data]]
  data[[argv$training_data]] <- NULL
}

debug(logger, paste0("Data Names: ", paste0(names(data), collapse = ", ")))

info(logger, "Training Lasso Model")
debug(logger, glue("LASSO Model column names: PrintGluedList(colnames(training))"))

# this gets weird
# Goal: build call using strings
# can not use eval_bare(parse_expr(str)) directly
# Writes call to temp file and then eval_bare(parse_expr(file))
# may be able to use base R, but I've tried too long, this works...
path <- tempfile()
if (argv$exclude) {
  debug(logger, "Excluding Vars")
  cat(
    "DoLassoModel(training, `",
    argv$training_axis,
    "`, exclude = ",
    argv$training_vars,
    ")",
    file = path,
    sep = ""
  )
} else {
  debug(logger, "Including Vars")
  cat("DoLassoModel(training, `",
    argv$training_axis,
    "`, include = ",
    argv$training_vars,
    ")",
    file = path,
    sep = ""
  )
}

debug(logger, glue("Lasso Model Command: {suppressWarnings(readLines(path))}"))

suppressWarnings({
  lasso_model <- rlang::parse_expr(file(path)) %>% rlang::eval_bare()
})

lasso_model_red <- coef(lasso_model)

lasso_model_red %>%
  as.data.frame() %>%
  mutate(rownames = rownames(lasso_model_red), .before = 1) %>%
  write_tsv(here(output_path, "lasso_model_coefs.tsv"))

lasso_model_red <- lasso_model_red[2:nrow(lasso_model_red), ]
lasso_model_red <- lasso_model_red[lasso_model_red != 0]

if (length(lasso_model_red) == 0) {
  error(logger, "The lasso model is empty!")
  source(here("lib/WriteInvocation.R"))
  WriteInvocation(argv, output_path = here(output_path, "invocation"))
  quit(status = 1)
}

times <- c(
  "Event Free Survival Time in Days",
  "Overall Survival Time in Days",
  "days_to_death",
  "OS_DAYS"
)
event_col <- c("First Event" = "relapse", "vital_status" = "Dead", "status" = 1)

debug(logger, "Doing Survival Analysis")
source(here("lib/batch_survival.R"))
suppressWarnings({
  results <- BatchSurvival(data, times, lasso_model, event_col)
})

pdf(here(plots_path, "survival.pdf"))
results$plots %>% walk(print)
graphics.off()

write_csv(results$stats, file = here(output_path, "survival.csv"))

source(here("lib/WriteInvocation.R"))
WriteInvocation(argv, output_path = here(output_path, "invocation"))
