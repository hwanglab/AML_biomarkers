#' Run a LASSO model
#' 
#' @param df a data.frame
#' @param outcome <data-masking> a column in df that has the continous outcome variable
#' @param exclusion <data-masking> a column (or vector of columns) in df to exclude from the lasso model
#' 
#' @return a \code{glmnet} object
#' 
#' @examples {
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
#' }
DoLassoModel <- function(df, outcome, exclude = NULL) {
  o <- rlang::enquo(outcome)
  e <- rlang::enquo(exclude)
  df <- as_tibble(df)
  
  if (rlang::quo_is_null(e)) {
    data <- dplyr::select(df, -!!o)
    response <- dplyr::pull(df, !!o)
    mat <- as.matrix(data)
  } else {
    cols <- tidyselect::eval_select(e, df)
    good_cols <- colnames(df)[colnames(df) %!in% names(cols)]
    selected <- rlang::set_names(df[-cols], good_cols)
    selected <- remove_missing(selected)
    data <- dplyr::select(selected, -!!o)
    response <- dplyr::pull(selected, !!o)
    mat <- as.matrix(data)
  }
  if (length(response) != nrow(mat)) rlang::abort("Number of observations does not match!")
  
  init_mod <- glmnet::cv.glmnet(mat, response)
  fin_mod <- glmnet::glmnet(mat, response, lambda = init_mod$lambda.min)
  return(fin_mod)
}

#' Calcualte scores using coefficents from a LASSO model
#' 
#' @param x a data.frame with colnames matching LASSO model
#' @param coef coeficents from a model
#' 
#' @return vector of scores
#' 
#' @examples {
#' # Not Run
#' # Using Base R
#' 
#' df$score <- UseLassoModelCoefs(df, coef(LASSO_model))
#' 
#' # Using dplyr
#' df %>%
#'    mutate(score = UseLASSOModelCoefs(cur_data(), coef(LASSO_model))
#' }
UseLASSOModelCoefs <- function(x, coef) {
  mod2 <- coef[-1, ] %>% as.data.frame()
  w_vec <- list(0)
  mod2 <- mod2[which(mod2[1] != 0), ,drop = FALSE]
  
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

#' Do A Kaplan-Meier Survival Analysis
#' 
#' @param data a dataframe (or coercible to a data.frame)
#' @param time,status,predictor <data-masking> varibles to compute survival on
#' @param method method for statstical testing
#' @param cox_ties method to use for ties in cox regression
#' @param group_by <data-masking> column to group survivial on
#' @param description A description of the data
#' @param lasso a lasso regression used during prediction
#' 
#' @return a list of stats and plots
#' 
#' @export
#' 
#' @description This function used dplyr style data-masking
#' 
DoSurvialAnalysis <- function(data, time, status, predictor,
                              method = "cox",
                              cox_ties = "efron",
                              group_by = NULL,
                              description = NULL,
                              lasso = NULL) {
  t <- rlang::eval_tidy(enquo(time), data = data)
  s <- rlang::eval_tidy(enquo(status), data = data)
  p <- rlang::eval_tidy(enquo(predictor), data = data)
  gb <- enquo(group_by)
  
  fit <- survminer::surv_fit(Surv(t, s) ~ p, data = data)
  
  surv_form <- survival::Surv(t, s) ~ p
  
  title <- paste0("Kaplan-Meier Curves for ",
                  rlang::as_name(enquo(status)),
                  " ~ ",
                  rlang::as_name(enquo(predictor)))
  xlab <- rlang::as_name(enquo(time))
  
  if(!is.null(description)) {
    subtitle1 <- paste0("Data: ", description)
  } else {
    subtitle1 <- NULL
  }
  
  if(!is.null(lasso)) {
    lasso <- coef(lasso) %>% as.matrix() 
    features <- lasso[lasso != 0, , drop = FALSE] %>% rownames()
    subtitle2 <- paste0("Model: ", paste(features[-1], collapse = ", "))
  } else {
    subtitle2 <- NULL
  }
  
  subtitle <- NULL
  if (!is.null(subtitle1) & !is.null(subtitle2)) {
    subtitle <- paste0(subtitle1, "\n", subtitle2)
  } else {
    if (!is.null(subtitle1)) subtitle <- subtitle1
    if (!is.null(subtitle2)) subtitle <- subtitle2
  }

  if (rlang::quo_is_null(gb)) {
    plot <- survminer::ggsurvplot(fit,
                       ylab = "Survival Probability",
                       pval = TRUE,
                       title = title,
                       xlab = xlab,
                       font.subtitle = 10,
                       subtitle = subtitle,
                       risk.table = "abs_pct")
  } else {
    g <- rlang::eval_tidy(gb, data = data)
    fit <- survminer::surv_fit(Surv(t, s) ~ g, data = data)
    plot <- survminer::ggsurvplot(fit,
                       ylab = "Survival Probability",
                       pval = FALSE,
                       title = paste0("Grouped ", title),
                       xlab = xlab,
                       font.subtitle = 8,
                       subtitle = subtitle,
                       risk.table = "abs_pct")
  }
  
  stat_res <- switch(EXPR = method,
                     "cox" = survival::coxph(surv_form, data, ties = cox_ties),
                     "chi-sq" = survival::survdiff(surv_form, data),
                     "none" = NULL,
                     rlang::abort("Unknown test: ", method))

  return(list(plot = plot, stat = stat_res))
} 

#' @param ensemble_col must be ensembl_gene_id_version (or others)
EnsemblToHGNC <- function(data, ensemble_col) {
  message("Accessing Biomart...")
  ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  ann <- biomaRt::getBM(attributes = c(ensemble_col,
                              "hgnc_symbol"),
                        mart = ensembl)
  
  message("Deduplicating counts...")
  data_named <- left_join(data, ann, by = ensemble_col) %>%
    remove_missing() %>%
    group_by(hgnc_symbol) %>%
    summarise(across(.cols = is_numeric, .fns = sum)) %>%
    filter(hgnc_symbol != "")
  return(data_named)
}

#' Transpose the data of a dataframe
#' 
#' @param df a data.frame
#' @param rowname_col <data-masking> column to call the observations, defaults to samples
#' @param colname_col <data-masking> column to turn into colnames, defaults to the first
#' 
#' @return a tibble with no rownames

transposeDF <- function(df, rowname_col = NULL, colname_col = NULL) {
  cn <- enquo(colname_col) %|q|% colnames(df)[[1]]
  rn <- enquo(rowname_col) %|q|% expr(samples)
  res <- df %>%
    select(-!!cn) %>%
    as.matrix() %>%
    t() %>%
    as_tibble() %>%
    magrittr::set_colnames(pull(df, !!cn)) %>%
    mutate(!!rn := colnames(df)[-1], .before = 1)
  return(res)
}

#' Test if Quosure is NULL
#' 
#' @concept `%||%` but for quosures
#' 
`%|q|%` <- function(x, y) {
  if (rlang::quo_is_null(x)) 
    y
  else x
}


#' Run CIBERSORTx
#' @param data
#' @param email
#' @param token
#' @param ref_sample
#' 

RunCIBERSORT <- function() {
  # test if docker is availible
  tryCatch(system("docker --help", ignore.stdout = TRUE),
           error = rlang::abort("Docker is not availible! CIBERSORTx requires docker. Please install docker."))
  
  data_file <- tempfile()
  write_tsv(data, file = data_file)
  
  # things to do:
  # write ref to tmp
  # enable doing ref
  # get token from enviorment var
  
  docker_cmd <- paste0("docker run -v ", tempdir(), ":/src/data ",
                       "-v ", tempdir(), ":/src/outdir cibersort/fractions ",
                       "--username ", email, " ",
                       "--token", token, " ",
                       "--refsample ", ref_sample, " ",
                       "--mixture ", data_file, " ",
                       "--fraction 0 --rmbatchSmode TRUE")
  
}

CalculatePValues <- function(x) {
  if (is.matrix(x$obs)) {
    otmp <- apply(x$obs, 1, sum)
    etmp <- apply(x$exp, 1, sum)
  }
  else {
    otmp <- x$obs
    etmp <- x$exp
  }
  
  df <- (sum(1 * (etmp > 0))) - 1
  pval <- pchisq(x$chisq, df, lower.tail = FALSE)
  return(pval)
}


  

