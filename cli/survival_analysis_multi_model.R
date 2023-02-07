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
  help = "how many values for alpha for the GLM should be tested? Alterntativly, a single value [0,1] to use as the alpha value.",
  default = 2
)
parser$add_argument(
  "--info",
  help = "a colon sepeated list including the name of the dataset, event, status, efs, etc",
  nargs = "+"
)
parser$add_argument(
  "--training",
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

renv::use_python(type = "conda")

source(here("cli/lib/utils.R"))
source(here("cli/lib/lasso.R"))
source(here("cli/lib/ml.R"))
source(here("cli/lib/stat_helpers.R"))

# source function from PR: https://github.com/kassambara/survminer/pull/582
source(here("cli/lib/ggforest.R"))

##### Prepare Directory --------------------------------------------------------
logger <- logger(threshold = argv$verbose)

if (argv$n_alpha > 0) {
  error(logger, "The value for n_alpha cannot be negative")
}

output_path <- PrepareOutDir(argv)
data <- LoadAnnotatedMatrix(output_path, logger)

output_path <- here(output_path, argv$test_id)
plots_path <- here(output_path, "plots")
info(logger, c("Using ", argv$test_id, " as the test id"))
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_path, recursive = TRUE, showWarnings = FALSE)

StopIfOutputDirNotExist(output_path)

if (argv$rerun_training) unlink(here(output_path, "saved_models"), recursive = TRUE)

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

debug(logger, "Preparing Data for training")
debug(logger, paste0("Data Names: ", paste0(names(data), collapse = ", ")))

training <- data[["TRAIN"]]
data[["TRAIN"]] <- NULL

### Parse Options for ML Training and Testing ----------------------------------
debug(logger, "Parsing Input options for ML")

params <- ParseNames(argv$info) %>%
  bind_rows() %>%
  mutate(
    cols = split(
      select(cur_data(), -contains("event"), -set),
      1:nrow(cur_data())
    ),
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
  inner_join(params, by = "set")

print(data_names)
### Prepare data for ML Methods ------------------------------------------------
debug(logger, "Preparing data")

training_params <- GetAndValidateTrainingParams(output_path, argv, training)
#print(training_params)
train_df <- do.call("CleanData", args = training_params)

training_params$df <- train_df
train_df <- do.call("PrepareDataForML", args = training_params)

non_cluster_names <- colnames(train_df) %>% str_subset("cluster", negate = TRUE)
cluster_ids <- colnames(train_df) %>%
  str_subset("cluster") %>%
  str_remove("cluster") %>%
  as.numeric() %>%
  range()
cluster_seq <- glue("clusters {cluster_ids[[1]]} - {cluster_ids[[2]]}")

collapsed_colnames <- glue_collapse(
  c(non_cluster_names, cluster_seq),
  sep = ", ",
  last = " and "
)
debug(logger, glue("The columns selected for training are: {collapsed_colnames}"))

## train_df <- CleanData(
##   training,
##   wbc = "wbc_at_diagnosis",
##   efs = "event_free_survival_time_in_days",
##   status = "First Event",
##   time_unit = "months"
## ) %>%
##   PrepareDataForML(wbc = "wbc_at_diagnosis")

# TODO: allow conversion of time to event from datasets.json
input_params <- select(data_names, data, cols)

data_names <- mutate(
  data_names,
  clean = pmap(
    input_params,
    ~ CleanData(..1, ..2, time_unit = "months")
  )
)

data_names <- mutate(data_names, testing = map(clean, PrepareDataForML, testing = FALSE))

if (length(argv$additional_testing_sets) != 0) {
  debug(logger, "Adding addional testing sets")
  add_set <- data_names %>%
    filter(value %in% argv$additional_testing_sets) %>%
    pull(testing)
  add_set <- append(add_set, list(train_df))
  train_df <- reduce(add_set, bind_rows)
}

##### LASSO Model --------------------------------------------------------------

info(logger, glue("Running Regularized GLM with {argv$n_alpha} value(s) of alpha"))
model_use <- "reg_glm"
if (CheckSavedModel(model_use)) {
  info(logger, glue("Existing model found. Using that instead..."))
  models <- LoadSavedModel(model_use)
} else {
  if (argv$n_alpha <= 1) {
    info(logger, "One value detected for alpha.")
    alphas <- argv$n_alpha
  } else {
    alphas <- seq(0, 1, length.out = as.integer(argv$n_alpha))
  }
  models <- DoGLMAlpha(
    train_df,
    "label",
    alphas,
    nfolds = as.numeric(argv$nfolds),
    k = 1000
  )
  saveRDS(models, GetSavedModelFilename(model_use))
  debug(logger, glue("{model_use} Model Training Done!"))
}

### SVM Model -----------------------------------------------------------------

info(logger, "Running SVM")
model_use <- "svm"
if (CheckSavedModel(model_use)) {
  info(logger, glue("Existing model found. Using that instead..."))
  models[[model_use]] <- LoadSavedModel(model_use)
  write_invoke <- FALSE
} else {
  models[[model_use]] <- svm(
    label ~ .,
    data = train_df,
    kernel = "polynomial",
    cost = 10
  )
  saveRDS(models[[model_use]], GetSavedModelFilename(model_use))
  debug(logger, glue("{model_use} Model Training Done!"))
}

### Random Forest Model --------------------------------------------------------

info(logger, "Running Random Forest")
model_use <- "rf"
write_invoke <- TRUE
if (CheckSavedModel(model_use)) {
  info(logger, glue("Existing model found. Using that instead..."))
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
    doBest = TRUE,
    plot = FALSE
  )
  saveRDS(models[[model_use]], GetSavedModelFilename(model_use))
  debug(logger, glue("{model_use} Model Training Done!"))
}

### KNN Regression Model -------------------------------------------------------

info(logger, "Running KNN Regression")
model_use <- "knn_regression"
if (CheckSavedModel(model_use)) {
  info(logger, glue("Existing model found. Using that instead..."))
  models[[model_use]] <- LoadSavedModel(model_use)
  write_invoke <- FALSE
} else {
  models[[model_use]] <- knnreg(label ~ ., data = train_df)
  saveRDS(models[[model_use]], GetSavedModelFilename(model_use))
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

models <- map(models, CheckLASSOValidity) %>% discard(is_null)

for (i in seq_along(models)) {
  debug(logger, glue("Running cox survival with {class(models[[i]])}"))
  palette <- c("#F9DC5C", "#011936")
  model_outs <- data_names %>%
    mutate(
      scores = map(clean, ~ CalcualteScoreFromModel(.x, models[[i]])),
      range = map(scores, ~ .x[["score"]]) %>% map_dbl(~ diff(range(.x)))
    ) %>%
    filter(range != 0, is.finite(range))

  if (nrow(model_outs) == 0) {
    info(logger, "Model created a range of 0 for all datasets")
    info(logger, "Skipping...")
  }
  model_outs <- model_outs %>%
    mutate(
      surv_data_bound = map2(surv_data, scores, bind_cols),
      fits = map(
        surv_data_bound,
        ~ survminer::surv_fit(Surv(time, code) ~ score_bin, data = .x)
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
          risk.table = "absolute",
          palette = palette
        )
      ),
      code_u = map_dbl(surv_data_bound, ~ length(unique(.x[["code"]]))),
      sco_bin_u = map_dbl(surv_data_bound, ~ length(unique(.x[["score_bin"]])))
    )
  notify_failures <- filter(model_outs, code_u <= 1 | sco_bin_u <= 1) %>%
    pull(value)
  model_outs <- filter(model_outs, code_u > 1, sco_bin_u > 1) %>%
    mutate(
      stat = map(
        surv_data_bound,
        ~ survival::coxph(Surv(time, code) ~ score_bin, data = .x)
      ),
      stat_out = map(stat, broom::tidy),
      liklihood_ratio = map_dbl(stat, AngrilyExtractLiklihoodRatio)
    )

  if (length(notify_failures) > 0) {
    excluded_sets <- glue_collapse(notify_failures, sep = ", ", last = ", ")
    warn(
      logger,
      glue(
        "The following datasets could not be dicotamized \\
        and are excluded: {excluded_sets}"
      )
    )
  }

  debug(logger, "Create New Directories for models")
  suppressWarnings({
    new_out_path <- glue("{output_path}/{names(models)[[i]]}")
    new_plot_path <- glue("{plots_path}/{names(models)[[i]]}")
  })

  dir.create(new_out_path, showWarnings = TRUE)
  dir.create(new_plot_path, showWarnings = TRUE)

  model_outs %>%
    select(where(is_character), liklihood_ratio, stat_out) %>%
    unnest_wider(stat_out) %>%
    write_tsv(glue("{new_out_path}/survival.tsv"))

  pdf(glue("{new_plot_path}/survival.pdf"))
  walk(pull(model_outs, plots), print)
  graphics.off()

  debug(logger, "Setting up TARGET Data...")
  surv_form_no_score <- Surv(time, code) ~ AR +
    `WT1 mutation` + `NPM mutation` + `CEBPA mutation`

  surv_form <- Surv(time, code) ~ score_bin + AR +
    `WT1 mutation` + `NPM mutation` + `CEBPA mutation`

  target_data <- filter(model_outs, value == "TARGET:FLT3")

  target_res <- mutate(target_data,
    forest_data = map2(data, surv_data_bound, SetupTARGETData)
  )

  fits1 <- map(target_res$forest_data, DoCox, surv_form)
  fits2 <- map(target_res$forest_data, DoCox, surv_form_no_score)

  fits1 <- transpose(fits1)
  fits2 <- transpose(fits2)

  target_res <- mutate(target_data,
    forest_data = map2(data, surv_data_bound, SetupTARGETData),
    stat1 = map(fits1[[1]], broom::tidy),
    stat2 = map(fits2[[1]], broom::tidy)
  )

  suppressMessages(
    target_res %>%
      select(stat1, where(is_character)) %>%
      unnest_auto(stat1) %>%
      unnest(cols = c(term, estimate, std.error, statistic, p.value)) %>%
      write_tsv(glue("{new_out_path}/cox_mutation_survival_with_score.tsv"))
  )

  pdf(glue("{new_plot_path}/cox_mutation_survival_with_score.pdf"))
  walk(fits1[[2]], print)
  graphics.off()
  suppressMessages(
    target_res %>%
      select(stat2, where(is_character)) %>%
      unnest_auto(stat2) %>%
      unnest(cols = c(term, estimate, std.error, statistic, p.value)) %>%
      write_tsv(glue("{new_out_path}/cox_mutation_survival_without_score.tsv"))
  )

  pdf(glue("{new_plot_path}/cox_mutation_survival_without_score.pdf"))
  walk(fits2[[2]], print)
  graphics.off()
}

### Write Invocation -----------------------------------------------------------
if (write_invoke) {
  source(here("lib/WriteInvocation.R"))
  WriteInvocation(argv, output_path = here(output_path, "invocation"))
}