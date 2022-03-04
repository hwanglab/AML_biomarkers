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
parser$add_argument(
  "--nfolds",
  "-F",
  help = "how many folds to use, 0 = LOO CV",
  default = 0
)
parser$add_argument(
  "--n-alpha",
  "-A",
  "-N",
  help = "how many values for alpha for the GLM should be tested? [currently only use 0]",
  default = 2
)
parser$add_argument(
  "--info",
  help = "a colon sepeated list including the name of the dataset, event, status, efs, etc",
  nargs = "+"
)
parser$add_argument(
  "--additional-testing-sets",
  "-T",
  help = "Additional annotated sets to mix into testing set",
  nargs = "+"
)
parser$add_argument(
  "--rerun-training",
  help = "Should training be redone?",
  action = "store_true"
)

argv <- parser$parse_args(c("-i", "flt3_cebpa_stem", "-I", "multi_test", "-E", "-v", "DEBUG", "-V", 'c(matches("^Year|^Age|^Overall"), where(is_character))', "-T", "AML02:FLT3", "--info", "AML02:efs=Event Free Survival Time in Days:wbc=WBC at Diagnosis:status=efs.evnt:event=1"))

# load more libraries
suppressPackageStartupMessages({
  library(here)
  library(log4r)
  library(survival)
  library(survminer)
  library(scorecard)
  library(glmnet)
  library(MASS)
  library(caret)
  library(janitor)
  library(randomForest)
  library(e1071)
  library(tensorflow)
  library(tfdatasets)
  library(keras)
  library(tidyverse)
  library(glue)
})

renv::use_python()

#### Make Functions ------------------------------------------------------------

#' Check the existence of a package
#
#' @param ... Package names
#' @param error If true, throw an error if the package doesn't exist
#
#' @return Invisibly returns boolean denoting if the package is installed
#'
#' @note Taken from Seurat
#
PackageCheck <- function(..., error = TRUE) {
  pkgs <- unlist(x = c(...), use.names = FALSE)
  package.installed <- vapply(
    X = pkgs,
    FUN = requireNamespace,
    FUN.VALUE = logical(length = 1L),
    quietly = TRUE
  )
  if (error && any(!package.installed)) {
    stop(
      "Cannot find the following packages: ",
      paste(pkgs[!package.installed], collapse = ", "),
      ". Please install"
    )
  }
  invisible(x = package.installed)
}

#' Not in operator
#'
#' @export
#' @rdname not_in
#' @name not_in
#'
#' @param lhs stuff
#' @param rhs stuff
#'

`%!in%` <- function(lhs, rhs) {
  match(lhs, rhs, nomatch = 0) == 0
}

#' Run a LASSO model
#'
#' @param df a data.frame
#' @param outcome <tidyselect> a column in df that has the continuous outcome variable
#' @param exclude <tidyselect> a column (or vector of columns) in df to exclude in the lasso model
#' @param include <tidyselect> a column (or vector of columns) in df to include in the lasso model
#' @param k number of LASSO regression iterations
#' @param nfolds number of folds for LASSO regression, 0 for Leave-one-out CV
#'
#' @return a \code{glmnet} object
#'
#' @importFrom rlang enquo set_names abort
#' @importFrom dplyr select pull
#' @importFrom tibble as_tibble

DoLassoModel <- function(df, outcome,
                         exclude = NULL,
                         include = NULL,
                         k = 1000,
                         nfolds = as.numeric(argv$nfolds)) {
  PackageCheck("glmnet")

  o <- rlang::enquo(outcome)
  e <- rlang::enquo(exclude)
  i <- rlang::enquo(include)
  df <- tibble::as_tibble(df)

  if (!rlang::quo_is_null(i) & !rlang::quo_is_null(e)) {
    rlang::abort("Both exclude and include cannot be provided")
  }

  NotNumeric <- function(data) {
    if (any(apply(data, 2, class) != "numeric")) {
      rlang::abort("Some columns are not numeric")
    }
  }

  # Subset the data based on include and exclude vars
  if (rlang::quo_is_null(e) & rlang::quo_is_null(i)) {
    NotNumeric(df)
    data <- dplyr::select(df, -!!o)
    response <- dplyr::pull(df, !!o)
    mat <- as.matrix(data)
  }
  if (!rlang::quo_is_null(i)) {
    # This is a way to do this
    data <- dplyr::select(df, !!i)
    cols <- tidyselect::eval_select(o, data, strict = FALSE) # gets column names
    data <- dplyr::select(data, -all_of(cols)) # removes those columns
    NotNumeric(data)
    response <- dplyr::pull(df, !!o)
    mat <- as.matrix(data)
  }
  if (!rlang::quo_is_null(e)) {
    cols <- tidyselect::eval_select(e, df)
    good_cols <- colnames(df)[colnames(df) %!in% names(cols)]
    selected <- rlang::set_names(df[-cols], good_cols)
    selected <- ggplot2::remove_missing(selected)
    NotNumeric(selected)
    data <- dplyr::select(selected, -!!o)
    response <- dplyr::pull(df, !!o)
    mat <- as.matrix(data)
  }

  if (length(response) != nrow(mat)) rlang::abort("Number of observations does not match!")
  if (nfolds == 0) nfolds <- length(response)
  fold_vector <- cv.glmnet(
    y = response,
    x = mat,
    nfolds = nfolds,
    keep = TRUE
  )$foldid

  LassoFunction <- function(alpha) {
    fit_res <- cv.glmnet(
      y = response,
      x = mat,
      foldid = fold_vector,
      alpha = alpha
    )
    fit_est <- as.numeric(coef(fit_res, s = "lambda.min"))
    return(fit_est)
  }

  AverageLassoModel <- function(x) {
    mask <- x != 0
    mask <- matrixStats::rowCounts(mask, value = TRUE)
    rows_select <- which(mask > 0.95 * k)
    B <- rowMeans(x[rows_select, ])
    return(B)
  }

  MultipleLassoFunction <- function(alpha) {
    info(logger, glue("    Running GLM with alpha = {alpha}"))
    map_dfc(1:k, ~ LassoFunction(alpha))
  }

  DetermineNSignificantFeatures <- function(x) {
    mask <- x != 0
    mask <- matrixStats::rowCounts(mask, value = TRUE)
    counts <- mask[mask != 0]
    return(length(counts))
  }

  alpha_cv <- tibble(alpha = seq(0, 1, length.out = as.numeric(argv$n_alpha)))

  # suppresses an annoying "New Names: " message
  suppressMessages({
    alpha_cv <- alpha_cv %>%
      mutate(
        model = map(alpha, MultipleLassoFunction),
        clean = map(model, AverageLassoModel),
        n     = map_int(model, DetermineNSignificantFeatures)
      )
  })

  cols_to_use <- colnames(data)
  mod_to_use <- filter(alpha_cv, alpha == 0) %>%
    pull(clean) %>%
    unlist(recursive = FALSE)
  names(mod_to_use) <- c("intercept", cols_to_use)

  return(mod_to_use)
}

#' Calcualte scores using coefficents from a LASSO model
#'
#' @param x a data.frame with colnames matching LASSO model
#' @param coef coeficents from a model
#'
#' @return vector of scores
#' @export
#'
#' @importFrom rlang warn
#' @importFrom purrr discard reduce
UseLASSOModelCoefs <- function(x, coef) {
  mod2 <- coef[-1, , drop = FALSE] # %>% as.matrix() %>% as.data.frame()
  w_vec <- list(0)
  mod2 <- mod2[which(mod2 != 0), , drop = FALSE]

  missing_features <- mod2[which(rownames(mod2) %!in% colnames(x)), , drop = FALSE] %>% rownames()
  n_features <- nrow(mod2)
  n_present <- n_features - length(missing_features)

  if (n_present == 0) rlang::abort(message = "There are no features matching the lasso model present!")

  if (length(missing_features) >= 1) {
    # build message
    mfeat <- "The following features are missing in the data: "
    feats <- paste(missing_features, collapse = ", ")
    pfeat <- paste0(
      "Only ",
      round(n_present / n_features * 100, digits = 1),
      "% features present."
    )
    rlang::warn(message = c(mfeat, feats, pfeat))
  }
  for (i in seq_along(rownames(mod2))) {
    col <- rownames(mod2)[i]
    vec <- x[[col]]
    w <- mod2[i, 1]
    w_vec[[i]] <- vec * w
  }
  w_vec <- purrr::discard(w_vec, is_empty)
  w_len <- length(w_vec)

  if (w_len == 1) {
    res <- w_vec[[1]]
  } else {
    res <- w_vec %>%
      purrr::reduce(rbind) %>%
      colSums()
  }
  return(res + coef[1, ])
}

SetupTARGETData <- function(df) {
  res <- df %>%
    mutate(
      event = if_else(`First Event` == "relapse", 1, 0),
      score_bin_median = if_else(score >= median(score), "LS-HIGH", "LS-LOW"),
      flt3_ratio = if_else(`FLT3/ITD allelic ratio` > 0.4, "AR+", "AR-"),
      MRD_end = if_else(`MRD % at end of course 1` > 5, "MRD+", "MRD-")
    )

  mutation_cols <- c(
    "WT1 mutation", "c-Kit Mutation Exon 8", "c-Kit Mutation Exon 17",
    "NPM mutation", "CEBPA mutation", "FLT3 PM", "MRD_end", "MLL"
  )

  res <- res %>%
    mutate(
      across(
        .cols = all_of(mutation_cols),
        .fns = ~ if_else(.x %in% c("not done", "unknown"), NA_character_, .x)
      )
    ) %>%
    mutate(
      `Cytogenetic Complexity` = if_else(
        `Cytogenetic Complexity` %in% c("na", "n/a"),
        NA_character_,
        `Cytogenetic Complexity`
      )
    )
  return(res)
}

GetSavedModelFilename <- function(type) {
  if (type == "nn_regression") {
    here(output_path, "saved_models", glue("{type}_model.h5"))
  } else {
    here(output_path, "saved_models", glue("{type}_model.rds"))
  }
}

CheckSavedModel <- function(type) {
  model_save <- GetSavedModelFilename(type)
  file.exists(model_save) || dir.exists(model_save)
}

LoadSavedModel <- function(type) {
  model_save <- GetSavedModelFilename(type)
  if (type == "nn_regression") {
    backend_k <- keras::backend()
    backend_k$clear_session()
    model <- load_model_hdf5(model_save)
  } else {
    model <- readRDS(model_save)
  }

  return(model)
}

##### Prepare Directory --------------------------------------------------------
logger <- logger(threshold = argv$verbose)

if (argv$dir == "") {
  output_path <- paste0("outs/", argv$id)
  plots_path <- paste0("plots/", argv$id)
} else {
  output_path <- paste0(parser$run_dir, "/outs/", argv$id)
  plots_path <- paste0(parser$run_dir, "/plots/", argv$id)
}

data_filename <- here(output_path, glue("cache/clinical_deconvoluted.rds"))

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
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
  dir.create(plots_path, recursive = TRUE, showWarnings = FALSE)
}

if (!dir.exists(here(output_path))) {
  fatal(logger, "Output directory does not exist")
  quit(status = 1)
}

if (argv$rerun) unlink(here(output_path, "saved_models"))

debug(logger, paste0("Importing Data from: ", data_filename))

suppressWarnings({
  data <- tryCatch(
    readRDS(data_filename),
    error = function(e) {
      error(logger, "Cannot find annotated deconvoluted samples.")
      quit(status = 1)
    }
  )
})

debug(logger, "Preparing Data for training")
debug(logger, paste0("Data Names: ", paste0(names(data), collapse = ", ")))

training <- data[["TRAIN"]]
data[["TRAIN"]] <- NULL

##### LASSO Model --------------------------------------------------------------
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
    ")\n",
    file = path,
    sep = ""
  )
}

write_invoke <- TRUE
if (CheckSavedModel("lasso")) {
  info(logger, "Existing LASSO model found. Importing Existing Model")
  lasso_model <- LoadSavedModel("lasso")
  write_invoke <- FALSE
} else {
  debug(logger, paste0("Lasso Model Command: ", suppressWarnings(readLines(path))))
  suppressWarnings({
    lasso_model <- rlang::parse_expr(file(path)) %>% rlang::eval_bare()
  })
  if (!dir.exists(here(output_path, "saved_models"))) {
    dir.create(here(output_path, "saved_models"))
  }
  saveRDS(lasso_model, GetSavedModelFilename("lasso"))
  debug(logger, "LASSO Model Training Done!")
}

lasso_model_red <- lasso_model %>%
  as.matrix() %>%
  as.data.frame()

lasso_model_red %>%
  rownames_to_column() %>%
  write_tsv(here(output_path, "lasso_model_coefs.tsv"))

lasso_model_red <- as.data.frame(lasso_model_red)
lasso_model_test <- lasso_model_red[2:nrow(lasso_model_red), , drop = FALSE]
lasso_model_test <- lasso_model_red[lasso_model_red != 0, , drop = FALSE]

if (length(lasso_model_test) == 0) {
  error(logger, "The lasso model is empty!")
  source(here("lib/WriteInvocation.R"))
  WriteInvocation(argv, output_path = here(output_path, "invocation"))
  quit(status = 1)
}

times <- c(
  "Event Free Survival Time in Days",
  "Overall Survival Time in Days",
  "days_to_death",
  "days_to_recurrance",
  "OS_DAYS"
)
event_col <- c("First Event" = "relapse", "vital_status" = "Dead", "status" = 1, "efs.evnt" = 1)

debug(logger, "Doing Survival Analysis")
source(here("lib/batch_survival.R"))
suppressWarnings({
  results <- BatchSurvival(data, times, lasso_model_red, event_col)
})

pdf(here(plots_path, "survival.pdf"))
results$plots %>% walk(print)
graphics.off()

write_csv(results$stats, file = here(output_path, "survival.csv"))

info(logger, "Mulitvariate analysis on TARGET")

data <- map2(data, results$scores, ~ mutate(.x, score = .y))
target_data <- data[c("TARGET:FLT3", "TARGET:NEG", "TARGET:CEBPA")]

debug(logger, "Setting up TARGET Data...")
target_data <- map(target_data, SetupTARGETData)
surv_form <- Surv(`Event Free Survival Time in Days`, event) ~ score_bin_median + flt3_ratio + MRD_end + `WT1 mutation`

palette <- c("#E7B800", "#2E9FDF")
fits <- coxph(surv_form, data = target_data[["TARGET:FLT3"]])

pdf(glue("{plots_path}/forest_plot.pdf"))
set <- "TARGET:FLT3"
fits <- coxph(surv_form, data = target_data[[set]])
ggforest(
  fits,
  data = target_data[[set]],
  main = glue("Hazard Ratios: {set}")
)

graphics.off()

### ML Functions  --------------------------------------------------------------
CleanData <- function(df, wbc = NULL, efs = NULL, status) {
  clean <- janitor::clean_names(df)

  wbc <- janitor::make_clean_names(wbc)
  efs <- janitor::make_clean_names(efs)
  status <- janitor::make_clean_names(status)

  namekey <- c("efs_days", "status")
  names(namekey) <- c(efs, status)

  clean <- plyr::rename(clean, namekey, warn_missing = TRUE) %>%
    select(starts_with("cluster"), efs_days, status)
  res <- clean %>%
    as_tibble(.name_repair = "minimal") %>%
    setNames(colnames(.)) %>%
    remove_missing()

  return(res)
}

PrepareDataForML <- function(df) {
  data <- dplyr::select(df, starts_with("cluster"))

  labs <- pull(df, efs_days)

  res <- data %>%
    as_tibble(.name_repair = "minimal") %>%
    setNames(colnames(.)) %>%
    mutate(label = labs)

  return(res)
}

PrepareDataForSurvival <- function(data, invert, event_use) {
  res <- dplyr::select(data, status, efs = efs_days)
  res$code <- if_else(res$status == event_use, 1, 0)
  if (invert) res$code <- if_else(res$status == event_use, 0, 1)
  return(res)
}

CalcualteScoreFromModel <- function(data, model) {
  score <- predict(model, data)
  data <- mutate(
    data,
    score = score,
    score_bin = if_else(score > median(score, na.rm = TRUE), "HIGH", "LOW")
  ) %>%
    dplyr::select(score, score_bin)
  return(data)
}

ReturnNamesFromString <- function() {
  if (length(argv$info) != 0) data_info <- c(data_info, argv$info)
  split <- str_split(data_info, ":")
  names <- map_chr(split, ~ .x[[1]])
  split1 <- map(split, ~ .x[-1])

  split2 <- map(split1, str_split, "=")
  names(split2) <- names
  map(
    split2,
    ~ map(.x, ~ set_names(.x[[2]], .x[[1]])) %>% unlist()
  ) %>%
    bind_rows() %>%
    mutate(set = names)
}

ReturnNamesFromJSON <- function() {
  if (!file.exists(argv$info)) {
    error(logger, glue("File does not exist: {argv$info}"))
  }
  jsonlite::fromJSON(argv$info) %>% map_dfr(as_tibble, .id = "set")
}

ParseNames <- function() {
  if (tools::file_ext(argv$info) == "json") {
    return(ReturnNamesFromJSON)
  } else {
    return(ReturnNamesFromString)
  }
}
### Parse Options for ML Training and Testing ----------------------------------
debug(logger, "Parsing Input options for ML")
data_info <- c(
  "TARGET:efs=Event Free Survival Time in Days:wbc=wbc_at_diagnosis:status=First Event:not_event=censored", # nolint
  "TCGA:efs=days_to_death:status=vital_status:event=Dead",
  "BeatAML:efs=OS_DAYS:status=status:event=1"
)

params <- ParseNames() %>%
  bind_rows() %>%
  mutate(
    invert = if_else(is.na(not_event), FALSE, TRUE),
    event_use = if_else(is.na(event), not_event, event)
  )

data_names <- names(data) %>%
  as_tibble() %>%
  mutate(
    set = str_split(value, ":"),
    set = map_chr(set, ~ .x[[1]]),
    data = data
  ) %>%
  full_join(params, by = "set")

### Prepare data for ML Methods ------------------------------------------------
debug(logger, "Preparing data")
train_df <- CleanData(
  training,
  wbc = "wbc_at_diagnosis",
  efs = "event_free_survival_time_in_days",
  status = "First Event"
) %>%
  PrepareDataForML()

input_params <- select(data_names, data, wbc, efs, status)
data_names <- mutate(data_names, clean = pmap(input_params, ~ CleanData(..1, ..2, ..3, ..4)))

data_names <- mutate(data_names, testing = map(clean, PrepareDataForML))

if (!is.na(argv$additional_testing_sets)) {
  debug(logger, "Adding addional testing sets")
  add_set <- filter(data_names, value %in% argv$additional_testing_sets) %>% pull(testing)
  add_set <- append(add_set, list(train_df))
  train_df <- reduce(add_set, bind_rows)
}

models <- list()
### SVM Model -----------------------------------------------------------------

info(logger, "Running SVM")
model_use <- "svm"
if (CheckSavedModel(model_use)) {
  info(logger, glue("Existing {model_use} model found. Importing Existing Model"))
  models[[model_use]] <- LoadSavedModel(model_use)
  write_invoke <- FALSE
} else {
  models[[model_use]] <- svm(label ~ ., data = train_df, kernel = "polynomial", cost = 10)
  saveRDS(models[[model_use]], GetSavedModelFilename(model_use))
  debug(logger, glue("{model_use} Model Training Done!"))
}

### Random Forest Model --------------------------------------------------------

info(logger, "Running Random Forests")
# varImpPlot(model)
# plot(model)
# sqrt(model$mse[which.min(model$mse)])

model_use <- "rf"
if (CheckSavedModel(model_use)) {
  info(logger, glue("Existing {model_use} model found. Importing Existing Model"))
  models[[model_use]] <- LoadSavedModel(model_use)
  write_invoke <- FALSE
} else {
  models[[model_use]] <- tuneRF(
    x = train_df %>% dplyr::select(-label),
    y = train_df$label,
    ntreeTry = 500,
    mtryStart = 13,
    stepFactor = 1.5,
    improve = 0.01,
    trace = FALSE,
    doBest = TRUE
  )
  saveRDS(models[[model_use]], GetSavedModelFilename(model_use))
  debug(logger, glue("{model_use} Model Training Done!"))
}

### KNN Regression Model -------------------------------------------------------

info(logger, "Running KNN Regression")
model_use <- "knn_regression"
if (CheckSavedModel(model_use)) {
  info(logger, glue("Existing {model_use} model found. Importing Existing Model"))
  models[[model_use]] <- LoadSavedModel(model_use)
  write_invoke <- FALSE
} else {
  models[[model_use]] <- knnreg(label ~ ., data = train_df)
  saveRDS(models[[model_use]], GetSavedModelFilename(model_use))
  debug(logger, glue("{model_use} Model Training Done!"))
}

### Neural Network Regression --------------------------------------------------

info(logger, "Running Neural Network Regression")

CallbackLogger <- callback_lambda(
  on_epoch_end = function(epoch, logs) {
    if (epoch %% 100 == 0 && epoch != 0) info(logger, glue("{epoch} epochs finished"))
  }
)
StopTrainingEarly <- callback_early_stopping(monitor = "val_loss", patience = 20)

BuildNN <- function(dataset, spec) {
  input <- layer_input_from_dataset(train_df %>% select(-label))

  output <- input %>%
    layer_dense_features(dense_features(spec_f)) %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_dense(units = 16, activation = "relu") %>%
    layer_dense(units = 1)

  model <- keras_model(input, output)

  model %>%
    compile(
      loss = "mse",
      optimizer = optimizer_nadam(learning_rate = 0.1),
      metrics = list("mean_absolute_error")
    )
  return(model)
}

spec_f <- feature_spec(train_df, !contains("label"), label) %>%
  step_numeric_column(all_numeric(), normalizer_fn = scaler_min_max()) %>%
  fit()

model_use <- "nn_regression"
if (CheckSavedModel(model_use)) {
  info(logger, glue("Existing {model_use} model found. Importing Existing Model"))
  models[[model_use]] <- LoadSavedModel(model_use)
  write_invoke <- FALSE
} else {
  spec_f <- feature_spec(train_df, !contains("label"), label) %>%
    step_numeric_column(all_numeric(), normalizer_fn = scaler_min_max()) %>%
    fit()
  models[[model_use]] <- BuildNN(train_df, spec_f)
  training_output <- models[[model_use]] %>%
    fit(
      x = train_df %>% select(-label),
      y = train_df$label,
      epochs = 100,
      validation_split = 0.2,
      verbose = 0,
      callbacks = list(CallbackLogger, StopTrainingEarly)
    )
  # save_model_hdf5(models[[model_use]], GetSavedModelFilename(model_use))
  debug(logger, glue("{model_use} Model Training Done!"))
}

## ML Survival Analysis --------------------------------------------------------

info(logger, "Doing Survival Analysis")

input_params <- dplyr::select(data_names, data = clean, invert, event_use)

data_names <- data_names %>%
  mutate(
    surv_data = pmap(input_params, PrepareDataForSurvival)
  )

for (i in seq_along(models)) {
  model_outs <- data_names %>%
    mutate(
      scores = map(clean, ~ CalcualteScoreFromModel(.x, models[[i]])),
      surv_data_bound = map2(surv_data, scores, bind_cols),
      fits = map(
        surv_data_bound,
        ~ survminer::surv_fit(Surv(efs, code) ~ score_bin, data = .x)
      ),
      plots = map2(
        fits,
        value,
        ~ survminer::ggsurvplot(.x,
          ylab = "Survival Probability",
          pval = FALSE,
          xlab = "Survival Time",
          font.subtitle = 8,
          subtitle = .y,
          risk.table = "abs_pct"
        )
      ),
      stat_out = map(
        surv_data_bound,
        ~ survival::coxph(Surv(efs, code) ~ score_bin, data = .x)
      ) %>%
        map(broom::tidy)
    )

  stat_res <- select(model_outs, stat_out, where(is_character)) %>% unnest_wider(stat_out)
  write_tsv(stat_res, glue("{output_path}/{names(models)[[i]]}_survival.tsv"))

  pdf(glue("{plots_path}/{names(models)[[i]]}_survival.pdf"))
  print(pull(model_outs, plots))
  graphics.off()
}

### Write Invocation -----------------------------------------------------------
if (write_invoke) {
  source(here("lib/WriteInvocation.R"))
  WriteInvocation(argv, output_path = here(output_path, "invocation"))
}
