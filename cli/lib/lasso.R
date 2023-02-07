

# Set custom classes for organization of objects for Lasso Modeling
setOldClass("cv.glmnet")
setOldClass("simpleError")
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
setClass("FailedModel", contains = c("character"))


#' Run a Lasso Model with predefined foldid and get coeffients
#' @param alpha alpha value to use
#' @param response response vector for prediction
#' @param mat feature matrix to predict with
#' @param fold_vector precomputed foldids
#'
#' @return GLM object
LassoFunction <- function(alpha, response, mat, fold_vector) {
  fit_res <- tryCatch(
    cv.glmnet(
      y = response,
      x = mat,
      foldid = fold_vector,
      alpha = alpha
    ),
  error = function(e) {
    ReturnFailedModel(e)
  }
  )

  if (is.FailedModel(fit_res)) return(fit_res)
  fit_est <- as.numeric(coef(fit_res, s = "lambda.min"))
  res <- new("GLM", glm = fit_res, coefs = fit_est)
  return(res)
}

#' Run k Lasso Models
#'
#' @inheritParams LassoFunction
#' @param k number of iterations for lasso model
MultipleLassoFunction <- function(alpha, response, mat, fold_vector, k) {
  purrr::map(1:k, ~ LassoFunction(alpha, response, mat, fold_vector)) %>%
    discard(is.FailedModel)
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
  response <- dplyr::pull(df, all_of(outcome))
  mat <- dplyr::select(df, -all_of(outcome)) %>% as.matrix()
  NotNumeric <- function(data) {
    if (any(apply(x, 2, class) != "numeric")) {
      rlang::abort("Some columns are not numeric")
    }
  }

  if (length(response) != nrow(mat)) rlang::abort("Number of observations does not match!")
  if (nfolds == 0) nfolds <- length(response)
  fold_vector <- tryCatch(
    cv.glmnet(
    y = response,
    x = mat,
    nfolds = nfolds,
    keep = TRUE
  )$foldid,
  error = function(e) {
    ReturnFailedModel(e)
  }
  )
  if (is.FailedModel(fold_vector)) {
    warning(logger, "GLM Model Training Failed, skipping...")
    return(fold_vector)
  }

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
DoGLMAlpha <- function(df, outcome, alphas, nfolds, k = 1000, save = TRUE) {
  res <- purrr::map(
    alphas,
    ~ DoLassoModel(df, outcome, alpha = .x, k = k, nfolds = nfolds)
  )
  names(res) <- glue::glue("LASSO_alpha_{alphas}")
  if (save) purrr::walk(res, SaveLassoAsTSV)
  return(res)
}

#' Save the Lasso coefficents to TSV file
#'
#' @param model a MultiLassoModel object
#'
#' @return invisible
SaveLassoAsTSV <- function(model) {
  if (is.FailedModel(model)) {
    return(invisible(NULL))
  }
  alpha <- model@alpha
  red <- model@average %>%
    `class<-`("data.frame") %>%
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
  if ("FailedModel" %in% class(model)) {
      model <- NULL
    }

  if (!("MultiLassoModel" %in% class(model))) {
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
      do.call(rbind, .) %>%
      as.data.frame() %>%
      colSums()
    names(res) <- NULL
  }
  return(res + object@average[1, ])
}

PrintGLM <- function(x, ...) {
  alpha <- x@alpha
  n_feats <- x@nfeats
  n_iter <- length(x@glms)
  p <- glue::glue("Multi-GLM Object with alpha of {alpha} and {n_feats} significant features ran {n_iter} times.")

  print(p)
}

#' Register method for predict
#' @method
.S3method("predict", "MultiLassoModel", UseLASSOModelCoefs)

#' Register method for print
#' @method
.S3method("print", "MultiLassoModel", PrintGLM)

setMethod(
  f = "show",
  signature = "MultiLassoModel",
  definition = function(object) {
    print(object)
  }
)

is.FailedModel <- function(x) inherits(x, "FailedModel")

ReturnFailedModel <- function(e) {
  message <- as.character(e)
  new("FailedModel", message)
}