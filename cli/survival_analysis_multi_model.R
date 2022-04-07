#!/usr/bin/env Rscript
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
  "--categorical-prediction",
  "-c",
  help = "Is the predictor value categorical?",
  action = "store_true"
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

source(here("cli/lib/utils.R"))
source(here("cli/lib/lasso.R"))
source(here("cli/lib/ml.R"))

##### Prepare Directory --------------------------------------------------------
logger <- logger(threshold = argv$verbose)

if (argv$dir == "") {
  output_path <- paste0("outs/", argv$id)
} else {
  output_path <- paste0(parser$run_dir, "/outs/", argv$id)
}

data_filename <- here(output_path, glue("cache/clinical_deconvoluted.rds"))

output_path <- here(output_path, argv$test_id)
plots_path <- here(output_path, "plots")
info(logger, c("Using ", argv$test_id, " as the test id"))
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_path, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(here(output_path))) {
  fatal(logger, "Output directory does not exist")
  quit(status = 1)
}

if (argv$rerun) unlink(here(output_path, "saved_models"))

if (dir.exists(here(output_path, "saved_modes"))) {
  if (argv$categorical) {
    if (!file.exists(here(output_path, "saved_models", "categorical"))) {
      error(logger, "Categorical Analysis requested, but ran in a Regression test directory!")
      quit(error = 1)
    }
  } else {
    if (!file.exists(here(output_path, "saved_models", "categorical"))) {
      error(logger, "Regression Analysis requested, but ran in a Categorical test directory!")
      quit(error = 1)
    }
  }
}

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

if (categorical) {
  models <- list()
} else {
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

if (argv$categorical) cat("", file = here(output_path, "saved_models", "categorical"))
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
