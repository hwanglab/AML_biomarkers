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

argv <- parser$parse_args()

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

# Set custom classes for organization of objects for Lasso Modeling
setOldClass("cv.glmnet")
setClass("AverageLassoModel", contains = "data.frame")
setClass("GLM", representation(glm = "cv.glmnet", coefs = "numeric"))
setClass(
  "MultiLassoModel",
  representation(
    glms = "list",
    average = "AverageLassoModel",
    nfeats = "numeric",
    alpha = "numeric"
  )
)

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
#' @param lhs stuff
#' @param rhs stuff
#'

`%!in%` <- function(lhs, rhs) {
  match(lhs, rhs, nomatch = 0) == 0
}

#' Run a Lasso Model with predefined foldid and get coeffients
#' @param alpha alpha value to use
#' @param response response vector for prediction
#' @param mat feature matrix to predict with
#' @param fold_vector precomputed foldids
#'
#' @return GLM object
LassoFunction <- function(alpha, response, mat, fold_vector) {
  fit_res <- cv.glmnet(
    y = response,
    x = mat,
    foldid = fold_vector,
    alpha = alpha
  )
  fit_est <- as.numeric(coef(fit_res, s = "lambda.min"))
  res <- new("GLM", glm = fit_res, coefs = fit_est)
  return(res)
}

#' Run k Lasso Models
#'
#' @inheritParams LassoFunction
#' @param k number of iterations for lasso model
MultipleLassoFunction <- function(alpha, response, mat, fold_vector, k) {
  purrr::map(1:k, ~ LassoFunction(alpha, response, mat, fold_vector))
}

#' Average Multiple iterations of a Lasso model
#'
#' @inheritParams LassoFunction
#' @inheritParams MultipleLassoFunction
#' @param objects A list of GLM object
#'
#' @return MultipleLassoModel object
AverageLassoModels <- function(objects, alpha, data, k) {
  x <- as.data.frame(purrr::map_dfc(objects, ~ .x@coefs))
  cols_to_use <- colnames(data)
  rownames(x) <- c("intercept", cols_to_use)
  mask <- x != 0
  mask <- matrixStats::rowCounts(mask, value = TRUE)
  rows_select <- which(mask > 0.95 * k)
  B <- rowMeans(x[rows_select, ]) %>% as.data.frame()
  B <- new("AverageLassoModel", B)
  nfeats <- DetermineNSignificantFeatures(objects)
  object <- new(
    "MultiLassoModel",
    glms = objects,
    average = B,
    nfeats = nfeats,
    alpha = alpha
  )
  return(object)
}

#' Find Number of features with non-zero weight
#' @param object A list of GLM object
#'
#' @return integer number of features in model
DetermineNSignificantFeatures <- function(objects) {
  x <- purrr::map_dfc(objects, ~ .x@coefs)
  mask <- x != 0
  mask <- matrixStats::rowCounts(mask, value = TRUE)
  counts <- mask[mask != 0]
  return(length(counts) - 1) # -1 for intercept
}

#' Run a LASSO model
#'
#' @param df a data.frame
#' @param outcome  a column in df that has the continuous outcome variable
#' @param k number of LASSO regression iterations
#' @param nfolds number of folds for LASSO regression, 0 for Leave-one-out CV
#'
#' @return MultiLassoModel
#'
#' @importFrom rlang enquo set_names abort
#' @importFrom dplyr select pull
#' @importFrom tibble as_tibble

DoLassoModel <- function(df, outcome,
                         alpha = 0,
                         k = 1000,
                         nfolds = 1) {
  PackageCheck("glmnet")
  df <- tibble::as_tibble(df)
  response <- dplyr::pull(df, outcome)
  mat <- dplyr::select(df, -outcome) %>% as.matrix()
  NotNumeric <- function(data) {
    if (any(apply(x, 2, class) != "numeric")) {
      rlang::abort("Some columns are not numeric")
    }
  }

  if (length(response) != nrow(mat)) rlang::abort("Number of observations does not match!")
  if (nfolds == 0) nfolds <- length(response)
  fold_vector <- cv.glmnet(
    y = response,
    x = mat,
    nfolds = nfolds,
    keep = TRUE
  )$foldid

  models <- MultipleLassoFunction(
    alpha = alpha,
    response = response,
    mat = mat,
    fold_vector = fold_vector,
    k = k
  )
  names(models) <- 1:k
  model <- AverageLassoModels(models, alpha = alpha, data = mat, k = k)
  return(model)
}

#' Do a GLM for each value of alpha
#'
#' @inheritParams DoLassoModel
#' @param alphas vector of alphas to use
#'
#' @return list of MultiLassoModel objects
DoGLMAlpha <- function(df, outcome, alphas, nfolds, k) {
  res <- purrr::map(
    alphas,
    ~ DoLassoModel(df, outcome, alpha = .x, k = k, nfolds = nfolds)
  )
  names(res) <- glue::glue("LASSO_alpha_{alphas}")
  purrr::walk2(res, alphas, SaveLassoAsTSV)
  return(res)
}

#' Save the Lasso coefficents to TSV file
#'
#' @param model a MultiLassoModel object
#'
#' @return invisible
SaveLassoAsTSV <- function(model) {
  alpha <- model@alpha
  red <- model@average %>%
    tibble::rownames_to_column() %>%
    rlang::set_names(c("feature", "avg_coef")) %>%
    write_tsv(here(glue("{output_path}/glm_coef_alpha_{alpha}.tsv")))
  return(invisible(NULL))
}

#' Check Validity of Lasso Model (2+ Features present)
#'
#' @param model An Object
#'
#' @return if model is not a MultiLassoModel object, return unchanged, \
#'         otherwise return object if the number of features is at least 2
CheckLASSOValidity <- function(model) {
  if (class(model) != "MultiLassoModel") {
    return(model)
  }
  if (model@nfeats == 0) model <- NULL
  return(model)
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
UseLASSOModelCoefs <- function(object, data = NULL) {
  if (is.null(data)) {
    rlang::abort(message = "Cannot predict with no data")
  }
  mod2 <- object@average[-1, , drop = FALSE]
  w_vec <- list()
  mod2 <- mod2[which(mod2 != 0), , drop = FALSE]

  missing_features <- mod2[which(rownames(mod2) %!in% colnames(data)), , drop = FALSE] %>% rownames()
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
    vec <- data[[col]]
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
  return(res + object@average[1, ])
}

#' Register method for predict
#' @method
.S3method("predict", "MultiLassoModel", UseLASSOModelCoefs)

#' Get the filename of a saved model
#' @param type a string used to identify the model
GetSavedModelFilename <- function(type) {
  if (type == "nn_regression") {
    here(output_path, "saved_models", glue("{type}_model.h5"))
  } else {
    here(output_path, "saved_models", glue("{type}_model.rds"))
  }
}

#' Check if model has been saved
#' @inheritParams GetSavedModelFilename
CheckSavedModel <- function(type) {
  model_save <- GetSavedModelFilename(type)
  file.exists(model_save) || dir.exists(model_save)
}

#' Load a saved model from disk
#' @inheritParams GetSavedModelFilename
LoadSavedModel <- function(type) {
  model_save <- GetSavedModelFilename(type)
  if (type == "nn_regression") {
    backend_k <- keras::backend()
    backend_k$clear_session()
    model <- keras::load_model_hdf5(model_save)
  } else {
    model <- readRDS(model_save)
  }
  return(model)
}

### ML Functions  --------------------------------------------------------------

#' Clean data to prepare for training and prediction
#'
#' @param df data to clean
#' @param wbc column with % WBCs as string
#' @param efs column with survival as string
#' @param status column with status as string
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

#' Subset data with good features
#' @param df a datafame
PrepareDataForML <- function(df) {
  data <- dplyr::select(df, starts_with("cluster"))

  labs <- pull(df, efs_days)

  res <- data %>%
    as_tibble(.name_repair = "minimal") %>%
    setNames(colnames(.)) %>%
    mutate(label = labs)

  return(res)
}

#' Recode status for prediction
#' @param data a dataframe
#' @param invert should the events be inverted (event_use are the good responses)
#' @param event_use what to call an event, should be 1 value in status column
PrepareDataForSurvival <- function(data, invert, event_use) {
  res <- dplyr::select(data, status, efs = efs_days)
  res$code <- if_else(res$status == event_use, 1, 0)
  if (invert) res$code <- if_else(res$status == event_use, 0, 1)
  return(res)
}

#' Calcuate a risk score from a model
#' @param data dataframe to use
#' @param model a model to use with a predict method
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

#' Parse information for survival
#'
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
    return(ReturnNamesFromJSON())
  } else {
    return(ReturnNamesFromString())
  }
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

if (!dir.exists(here(output_path, "saved_models"))) {
  dir.create(here(output_path, "saved_models"))
}

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

##### LASSO Model --------------------------------------------------------------
info(logger, "Training Lasso Model")

info(logger, "Running LASSO")
model_use <- "lasso"
if (CheckSavedModel(model_use)) {
  info(logger, glue("Existing {model_use} model found. Importing Existing Model"))
  models <- LoadSavedModel(model_use)
} else {
  alphas <- seq(0, 1, length.out = argv$n_alpha)
  models <- DoGLMAlpha(
    train_df,
    "label",
    alphas,
    nfolds = as.numeric(argv$nfolds)
  )
  saveRDS(models, GetSavedModelFilename(model_use))
  debug(logger, glue("{model_use} Model Training Done!"))
}

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

## Survival Analysis --------------------------------------------------------

info(logger, "Doing Survival Analysis")
input_params <- dplyr::select(data_names, data = clean, invert, event_use)

data_names <- data_names %>%
  mutate(
    surv_data = pmap(input_params, PrepareDataForSurvival)
  )

models <- map(models, CheckLASSOValidity)

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
      code_u = map_dbl(surv_data_bound, ~ length(unique(.x[["code"]]))),
      sco_bin_u = map_dbl(surv_data_bound, ~ length(unique(.x[["score_bin"]])))
    )
  model_outs <- filter(model_outs, code_u > 1, sco_bin_u > 1) %>%
    mutate(
      stat_out = map(
        surv_data_bound,
        ~ survival::coxph(Surv(efs, code) ~ score_bin, data = .x)
      ) %>%
        map(broom::tidy)
    )

  stat_res <- select(model_outs, stat_out, where(is_character)) %>%
    unnest_wider(stat_out)
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
