#!/usr/bin/env Rscript
source("renv/activate.R")
library(argparse)

# parse args
parser <- ArgumentParser("Do Survival Analysis")
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
  "--test-id",
  "-I",
  help = "ID to use for outputs. Must be unique.",
  default = "incremental"
)
parser$add_argument(
  "--verbose",
  "-v",
  help = "should messages be printed? One of: DEBUG, INFO, WARN, ERROR",
  default = "INFO"
)
parser$add_argument(
  "--sig-clusters",
  "-s",
  help = "should significant clusters be filtered on?",
  action = "store_true"
)
parser$add_argument(
  "--split",
  "-S",
  help = "should a saved file be read in to split patients?",
  action = "store_true"
)
parser$add_argument(
  "--training-data",
  "-t",
  help = "what should be used to train on? (FLT3, CEBPA, NEG, BeatAML, TCGA)",
  default = "FLT3"
)
parser$add_argument(
  "--training-split",
  "-T",
  help = "should the training data be split? 0 = no split",
  default = 0.7
)
parser$add_argument(
  "--training-axis",
  "-a",
  help = "What time axis to use for training",
  default = "Event Free Survival Time in Days"
)
parser$add_argument(
  "--training-vars",
  "-V",
  help = "What vars to include during Training. Should be a quoted R expression"
)
parser$add_argument(
  "--exclude",
  "-E",
  help = "Should training-vars be excluded instead?",
  action = "store_true"
)

argv <- parser$parse_args()

# load more libraries
suppressPackageStartupMessages({
  library(here)
  library(log4r)
  library(survival)
  library(survminer)
  library(scorecard)
  library(glmnet)
  suppressWarnings(library(rsinglecell))
  library(tidyverse)
})

logger <- logger(threshold = argv$verbose)

if (argv$dir == "") {
  output_path <- paste0("outs/", argv$id)
  plots_path <- paste0("plots/", argv$id)
} else {
  output_path <- paste0(parser$run_dir, "/outs/", argv$id)
  plots_path <- paste0(parser$run_dir, "/plots/", argv$id)
}

data_filename <- here(output_path, "cache/clinical_deconvoluted.rds")

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

debug(logger, paste0("Importing Data from: ", data_filename))

data <- tryCatch(
  readRDS(data_filename),
  error = function(e) {
    error(logger, "Cannot find annotated deconvoluted samples.")
    quit(status = 1)
  }
)

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
  saveRDS(
    map(split, ~ .x[[1]]),
    file = here("outs", argv$id, "cache", paste0(argv$training_data, "_rownames.rds"))
  )
}

debug(logger, "Preparing Data for training")

if (argv$training_split != 0) {
  training <- split$train
  data[[argv$training_data]] <- NULL
  data[[paste0(argv$training_data, ": Training")]] <- split$training
  data[[paste0(argv$training_data, ": Test")]] <- split$test
} else {
  training <- data[[argv$training_data]]
  data[[argv$training_data]] <- NULL
}

debug(logger, paste0("Data Names: ", paste0(names(data), collapse = ", ")))

if (argv$sig_clusters) {
  sig_cluster_file <- here(output_path, "cache/survival_models.rds")
  if (!file.exists(sig_cluster_file)) {
    error(logger, "Significant clusters have not been found")
  }
  debug(logger, "Filtering by significant clusters")
  survival_models <- readRDS(sig_cluster_file)
  sig_clusters <- survival_models %>%
    filter(chi_sig == "chi_sig") %>%
    pull(cluster) %>%
    paste0("cluster", .)
  training <- select(training, all_of(sig_clusters), matches("^[^cluster]"))
}

info(logger, "Training Lasso Model")

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

debug(logger, paste0("Lasso Model Command: ", suppressWarnings(readLines(path))))

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
