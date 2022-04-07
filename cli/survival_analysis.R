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
  "--sig-clusters",
  "-s",
  help = "should significant clusters be filtered on?",
  action = "store_true"
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
argv <- parser$parse_args()

# load more libraries
suppressPackageStartupMessages({
  library(here)
  library(log4r)
  library(survival)
  library(survminer)
  library(scorecard)
  library(glmnet)
  library(tidyverse)
  library(glue)
})


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
      paste(pkgs[!package.installed], collapse = ', '),
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
#' 
#' @export
#' 
#' @examples {
#' \dontrun{
#' DoLassoModel(target_split$train, `Event Free Survival Time in Days`)
#' 
#' # one exclusion column
#' target_model <- DoLassoModel(target_split$train,
#'                              outcome = `Event Free Survival Time in Days`,
#'                              exclude = `First Event`)
#' # multiple exclusion columns
#' DoLassoModel(target_split$train, 
#'              outcome = `Event Free Survival Time in Days`,
#'              exclude = c(`First Event`, AAR2, `FLT3/ITD allelic ratio`))
#' DoLassoModel(lasso_data,
#'              outcome = `Event Free Survival Time in Days`,
#'              exclude = c(`P-value`:RMSE, where(is_character)))
#'              }
#' }
DoLassoModel <- function(
    df, outcome,
    exclude = NULL,
    include = NULL,
    k = 1000,
    nfolds = as.numeric(argv$nfolds)
  ) {
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
#' 
#' @examples {
#' \dontrun{
#' # Using Base R
#' 
#' df$score <- UseLassoModelCoefs(df, coef(LASSO_model))
#' 
#' # Using dplyr
#' df %>%
#'    mutate(score = UseLASSOModelCoefs(cur_data(), coef(LASSO_model))
#' }
#' }
UseLASSOModelCoefs <- function(x, coef) {
  mod2 <- coef[-1, , drop = FALSE] #%>% as.matrix() %>% as.data.frame()
  debug(logger, "        Cleaning up model")
  w_vec <- list(0)
  mod2 <- mod2[which(mod2 != 0), ,drop = FALSE]
  
  missing_features <- mod2[which(rownames(mod2) %!in% colnames(x)), , drop = FALSE] %>% rownames()
  n_features <- nrow(mod2)
  n_present <- n_features - length(missing_features)
  
  if (n_present == 0) rlang::abort(message = "There are no features matching the lasso model present!")
  
  if (length(missing_features) >= 1)  {
    # build message
    mfeat <- "The following features are missing in the data: "
    feats <- paste(missing_features, collapse = ", ")
    pfeat <- paste0("Only ",
                    round(n_present / n_features * 100, digits = 1),
                    "% features present.")
    rlang::warn(message = c(mfeat, feats, pfeat))
  }
  debug(logger, "        Calculating weight vector")
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
    res <- w_vec %>% purrr::reduce(rbind) %>% colSums()
  }
  return(res + coef[1, ])
}
#' Calculate P values from a survival analysis
#' 
#' @param x result from either a chi square or cox regression
#' 
#' @return p value
#' @rdname pval-calc
#' @export

##### End Functions ------------------------------------------------------------
logger <- logger(threshold = argv$verbose)

if (argv$dir == "") {
  output_path <- paste0("outs/", argv$id)
} else {
  output_path <- paste0(parser$run_dir, "/outs/", argv$id)
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
  dir.create(output_path, showWarnings = FALSE)
} else {
  output_path <- here(output_path, argv$test_id)
  info(logger, c("Using ", argv$test_id, " as the test id"))
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
}

if (!dir.exists(here(output_path))) {
  fatal(logger, "Output directory does not exist")
  quit(status = 1)
}

plots_path <- here(output_path, "plots")
dir.create(plots_path, recursive = TRUE, showWarnings = FALSE)

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
    ")\n",
    file = path,
    sep = ""
  )
}

debug(logger, paste0("Lasso Model Command: ", suppressWarnings(readLines(path))))

lasso_model_rds <- here(output_path, glue("lasso_model.rds"))

write_invoke <- TRUE
if (!file.exists(lasso_model_rds)) {
  suppressWarnings({
    lasso_model <- rlang::parse_expr(file(path)) %>% rlang::eval_bare()
  })
  saveRDS(lasso_model, lasso_model_rds)
} else {
  info(logger, "Existing LASSO model found. Importing Existing Model")
  lasso_model <- readRDS(lasso_model_rds)
  write_invoke <- FALSE
}

debug(logger, "LASSO Model Training Done!")

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
event_col <- c("First Event" = "relapse", "vital_status" = "Dead", "status" = 1)

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
target_data <- data[c("FLT3", "NEG", "CEBPA")]

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

debug(logger, "Setting up TARGET Data...")
target_data <- map(target_data, SetupTARGETData)

surv_form <- Surv(`Event Free Survival Time in Days`, event) ~ 
  score_bin_median + flt3_ratio + `WT1 mutation` + `c-Kit Mutation Exon 8` + 
  `c-Kit Mutation Exon 17` + `NPM mutation` + `CEBPA mutation` + `MLL`

palette <- c("#E7B800", "#2E9FDF")
fits <- coxph(surv_form, data = target_data$`FLT3`)

pdf(glue("{plots_path}/forest_plot.pdf"))
set <- "FLT3"
fits <- coxph(surv_form, data = target_data[[set]])
ggforest(
  fits,
  data = target_data[[set]], 
  main = glue("Hazard Ratios: {set}")
  )

graphics.off()

if (write_invoke) {
  source(here("lib/WriteInvocation.R"))
  WriteInvocation(argv, output_path = here(output_path, "invocation"))
}
