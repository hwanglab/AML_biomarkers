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
  "--training-axis",
  "-a",
  help = "What time axis to use for training",
  default = "Event Free Survival Time in Days"
)
parser$add_argument(
  "--gene-set",
  "-G",
  help = "gene set name(s) to pull genes for",
  nargs = "+"
)

argv <- parser$parse_args()

# load more libraries
suppressPackageStartupMessages({
  library(clusterProfiler)
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

source("cli/lib/lasso_funs.R")

logger <- logger(threshold = argv$verbose)

source(here("cli/lib/utils.R"))
output_path <- PrepareOutDir(argv)
StopIfOutputDirNotExist(output_path)

bc_ext <- "no_bc"
if (argv$batch_correct) bc_ext <- argv$batch_correct_method

plots_path <- here(output_path, "plots")

# Plan
# get pathway genes
# prepare LASSO using those genes each individually
# assess univariate survival
# assess multivariate survival

ListGeneSets <- function(df) {
  gene_signatures <- dplyr::select(df, gene_symbol, gs_name)
  list <- dplyr::nest_by(gene_signatures, gs_name)
  names <- list[[1]]
  names(list[[2]]) <- names
  data <- purrr::map(list[[2]], c, use.names = TRUE)
  data <- unlist(data, recursive = FALSE)
  names(data) <- names
  return(data)
}

gene_sets <- msigdbr::msigdbr() %>% ListGeneSets()

gene_sets_selected <- gene_sets[argv$gene_set]

dataset_names <- c("target", "tcga", "beat_aml")

datasets <- map(
  dataset_names,
  ~ read_tsv(
    glue("cibersort_in/{.x}_data.txt"),
    col_types = cols()
  )
)

source("cli/lib/prepare_clin_info.R")

names(datasets) <- dataset_names
datasets <- map(datasets, rsinglecell::transposeDF)

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

debug(
  logger,
  glue(
    "TARGET has {nrow(annotated_datasets[['target']])} rows \\
    and {ncol(annotated_datasets[['target']])} cols"
  )
)

data_filename <- here(output_path, glue("cache/{bc_ext}_clinical_deconvoluted.rds"))

debug(logger, paste0("Importing Data from: ", data_filename))

cibersort_data <- tryCatch(
  readRDS(data_filename),
  error = function(e) {
    error(logger, "Cannot find annotated deconvoluted samples.")
    quit(status = 1)
  }
)

left_join_pos <- possibly(left_join, otherwise = NULL)
annotated_datasets <- map(annotated_datasets, ~ mutate(.x, across(where(is_character), str_to_lower)))
data_crossed <- purrr::transpose(cross2(cibersort_data, annotated_datasets))
data_crossed_names <- purrr::transpose(
  cross2(
    names(cibersort_data),
    names(annotated_datasets)
  )
)

cibersort_data_ann <- map2(data_crossed[[1]], data_crossed[[2]], left_join_pos)

names(cibersort_data_ann) <- map2(
  data_crossed_names[[1]],
  data_crossed_names[[2]],
  ~ glue("{.x} {.y}")
)
cibersort_data_ann <- discard(cibersort_data_ann, is_null)

# why.....
cibersort_data_ann[["BeatAML target"]] <- NULL
names(cibersort_data_ann) <- names(cibersort_data_ann) %>%
  str_extract(boundary("word"))

CreatePathwayModel <- function(data, genes, response_var = "Event Free Survival Time in Days") {
  pred <- select(data, any_of(genes))
  resp <- data[[response_var]]
  cv_mod <- cv.glmnet(as.matrix(pred), resp)
  min_cv <- cv_mod$lambda.min
  fin_mod <- glmnet(as.matrix(pred), resp, lambda = min_cv)
  return(coef(fin_mod))
}

# lasso_coefs is a list of lasso models from CreatePathwayModel
ReturnModelScores <- function(data, lasso_coefs) {
  PossiblyUseLMC <- possibly(UseLASSOModelCoefs, otherwise = 0)
  scores <- map_dfc(lasso_coefs, ~ PossiblyUseLMC(data, .x))
  colnames(scores) <- argv$gene_set
  return(cbind(data, scores))
}
training <- cibersort_data_ann$TRAIN
cibersort_data_ann$TRAIN <- NULL
models <- map(gene_sets_selected, ~ CreatePathwayModel(training, .x))
data_model_res <- map(cibersort_data_ann, ReturnModelScores, lasso_coefs = models)

training_model_res <- ReturnModelScores(training, lasso_coefs = models)

## Lets just remove TCGA and BeatAML for now...

data_model_res[["BeatAML"]] <- NULL
data_model_res[["TCGA"]] <- NULL

times <- c("Event Free Survival Time in Days", "Overall Survival Time in Days", "days_to_death", "OS_DAYS")

event_col <- c("First Event" = "relapse", "vital_status" = "Dead", "status" = 1)

library(rpart)



form <- reformulate(argv$gene_set[[1]], "`Event Free Survival Time in Days`")
test_rpart <- rpart(form, data = training_model_res, method = "class")

min_cp <- test_rpart$cptable[which.min(test_rpart$cptable[, "xerror"]), "CP"]

test_rpart_prune <- prune(test_rpart, cp = min_cp)

cutoffs <- test_rpart_prune$functions$print(
  test_rpart_prune$frame$yval2,
  attr(test_rpart_prune, "ylevels"),
  getOption("digits")
) %>%
  str_split(" ") %>%
  unlist()
cutoff <- cutoffs[1]

data_model_res$FLT3 <- mutate(
  data_model_res$FLT3,
  score_bin_rpart = if_else(ZZZ3_TARGET_GENES > cutoff, "HIGH", "LOW"),
  score_bin_median = if_else(ZZZ3_TARGET_GENES > median(ZZZ3_TARGET_GENES), "HIGH", "LOW"),
  event = if_else(`First Event` == "relapse", TRUE, FALSE)
)


form_rpart <- Surv(`Event Free Survival Time in Days`, event) ~ score_bin_median
coxph(form_rpart, data = data_model_res$FLT3)
palette <- c("#E7B800", "#2E9FDF")

cutoffs <- map(argv$gene_set, .f = function(pathway) {
  rpart_form <- reformulate(pathway, "`Event Free Survival Time in Days`")
  rpart_res <- rpart(form, data = training_model_res, method = "class")
  min_cp <- rpart_res$cptable[which.min(rpart_res$cptable[, "xerror"]), "CP"]

  prune_rpart <- prune(rpart_res, cp = min_cp)
  cutoffs <- prune_rpart$functions$print(
    prune_rpart$frame$yval2,
    attr(prune_rpart, "ylevels"),
    getOption("digits")
  ) %>%
    str_split(" ") %>%
    unlist()

  cutoff <- cutoffs[1]
  return(cutoff)
})

names(cutoffs) <- argv$gene_set


palette <- c("#E7B800", "#2E9FDF")

for (i in seq_along(cutoffs)) {
  pathway <- names(cutoffs)[[i]]
  cutoff <- cutoffs[[i]]
  for (dataset in data_model_res) {
    DifferentiateSamples <- function(x, value) {
      return(if_else(x > value, "HIGH", "LOW"))
    }
    work <- mutate(
      dataset,
      rpart  = DifferentiateSamples(.data[[pathway]], cutoff),
      median = DifferentiateSamples(.data[[pathway]], median(.data[[pathway]])),
      event = if_else(`First Event` == "relapse", TRUE, FALSE)
    )
    # Use RPART
    form <- Surv(`Event Free Survival Time in Days`, event) ~ rpart
    c1 <- coxph(form, data = work)
    fit <- surv_fit(form, data = work)
    p1 <- ggsurvplot(fit, data = work, palette = palette, conf.inv = TRUE)

    # use median
    form <- Surv(`Event Free Survival Time in Days`, event) ~ median
    c2 <- coxph(form, data = work)
    fit <- surv_fit(form, data = work)
    p2 <- ggsurvplot(fit, data = work, palette = palette, conf.inv = TRUE)
  }
}
