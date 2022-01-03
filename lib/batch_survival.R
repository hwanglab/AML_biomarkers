#' Do a batch survival analysis
#' 
#' @param data named list of datasets to test on
#' @param time columns where survival times are stored
#' @param lasso_model a lasso model
#' @param event_col a named vector where the names are columns and the values are strings of desired outcomes
#' 
#' @return a list with plots and data.frame of statitics
#' 
#' @examples 
#' /dontrun{
#' # setup times in the data
#' times <- c("Event Free Survival Time in Days", "Overall Survival Time in Days", "days_to_death", "OS_DAYS")
#' 
#' # set up what the event colums are and what the events are labelled as
#' event_col <- c("First Event" = "relapse", "vital_status" = "Dead", "status" = 1)
#' 
#' # run the survival analysis
#' results <- BatchSurvival(data, times, lasso_model, event_col)
#' }
BatchSurvival <- function(data, time, lasso_model, event_col) {

    events_ <- rlang::syms(names(event_col))
    times <- c("Event Free Survival Time in Days", "Overall Survival Time in Days", "days_to_death", "OS_DAYS")

    stats_final <- list(0)
    plots_final <- list(0)
    score_final <- list(0)
    for (dataset in seq_along(data)) {
        debug(logger, paste0("Batch Survival on: ", names(data)[[dataset]]))
        x <- as_tibble(data[[dataset]])
        avail_cols <- events_[names(event_col) %in% colnames(data[[dataset]])]
        avail_events <- event_col[names(event_col) %in% colnames(data[[dataset]])]

        stats <- list(0)
        plots <- list(0)
        
        for (st in seq_along(avail_cols)) {
            debug(logger, paste0("   Using ", avail_cols[[st]], " on ", names(data)[[dataset]]))
            x[["status"]] <- if_else(x[[avail_cols[[st]]]] == avail_events[[st]], 1, 0)
            debug(logger, "   Using Lasso coefs to calculate module score")
            x[["score"]] <- UseLASSOModelCoefs(x, coef(lasso_model))
            debug(logger, "   Binning scores on per dataset median")
            x[["score_bin"]] <- if_else(x$score >= median(x$score, na.rm = TRUE), "High", "Low")
            res <- list(0)

            for (time_axis in seq_along(times)) {
              debug(logger, paste0("      Using ", times[[time_axis]]))
                pSuv <- purrr::possibly(DoSurvivalAnalysis, otherwise = NA)
                desc <- paste0(names(data)[[dataset]], ": ", avail_cols[[st]])
                res[[time_axis]] <- pSuv(x,
                    !!rlang::sym(times[[time_axis]]),
                    status,
                    score,
                    group_by = score_bin,
                    description = desc,
                    lasso = lasso_model)
            }
          debug(logger, "Done with times!")
            names(res) <- times
            res <- purrr::discard(res, ~ all(is.na(.x)))
          debug(logger, paste0("Lenght of res is: ", length(res)))
            if (length(res) > 0) {
                stats_unclean <- map(res, ~ .x[["stat"]])
                plots[[st]] <- map(res, ~ .x[["plot"]])
                stats[[st]] <- stats_unclean %>%
                    map(broom::tidy, conf.int = TRUE) %>%
                    reduce(bind_rows) %>%
                    mutate(meta = names(res),
                        data = names(data)[[dataset]],
                        col = rep_len(rlang::as_name(avail_cols[[st]]), length(stats_unclean)))
            }
        }
        debug(logger, "Done with cols!")
        plots_final[[dataset]] <- plots
        stats_final[[dataset]] <- reduce(stats, bind_rows)
        score_final[[dataset]] <- x[["score"]]
    }
    debug(logger, "Done with data!")
    res1 <- stats_final %>% keep(is.data.frame) %>% reduce(bind_rows)
    debug(logger, "res1 created")
    res2 <- plots_final %>% unlist(recursive = FALSE) %>% unlist(recursive = FALSE)
    debug(logger, "res2 created")
    return(list(stats = res1, plots = res2, scores = score_final))
}

#' Do A Kaplan-Meier Survival Analysis
#' 
#' @param data a dataframe (or coercible to a data.frame)
#' @param time,status,predictor <data-masking> variables to compute survival on
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
#' @importFrom rlang eval_tidy enquo as_name
#' 
DoSurvivalAnalysis <- function(data, time, status, predictor,
                              method = "cox",
                              cox_ties = "efron",
                              group_by = NULL,
                              description = NULL,
                              lasso = NULL) {
  
  t <- rlang::eval_tidy(rlang::enquo(time), data = data)
  s <- rlang::eval_tidy(rlang::enquo(status), data = data)
  p <- rlang::eval_tidy(rlang::enquo(predictor), data = data)
  gb <- rlang::enquo(group_by)
  
  fit <- survminer::surv_fit(Surv(t, s) ~ p, data = data)
  
  surv_form <- survival::Surv(t, s) ~ p
  
  stat_res <- switch(EXPR = method,
                     "cox" = survival::coxph(surv_form, data, ties = cox_ties),
                     "chi-sq" = survival::survdiff(surv_form, data),
                     "none" = NULL,
                     rlang::abort("Unknown test: ", method))
  
  p_val <- switch(EXPR = method,
                  "cox" = CalculatePValues.cox(stat_res),
                  "chi-sq" = CalculatePValues.chi(stat_res),
                  "none" = NULL)
  
  title <- paste0("Kaplan-Meier Curves for ",
                  rlang::as_name(enquo(status)),
                  " ~ ",
                  rlang::as_name(enquo(predictor)))
  xlab <- rlang::as_name(enquo(time))
  
  if (!is.null(description)) {
    subtitle1 <- paste0("Data: ", description)
  } else {
    subtitle1 <- NULL
  }
  
  if (!is.null(lasso)) {
    lasso <- coef(lasso)[-1, , drop = FALSE] %>% as.matrix() %>% as.data.frame()
    features <- lasso[lasso != 0, , drop = FALSE] %>%
      dplyr::mutate(x = abs(s0)) %>%
      dplyr::arrange(dplyr::desc(x)) %>% rownames()
    if (length(features) >= 8) features <- c(features[1:8], "...")
    subtitle2 <- paste0("Model: ", paste(features, collapse = ", "))
    
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
  
  if (method != "none") subtitle <- paste0(subtitle, "\np = ", p_val)
  
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
  return(list(plot = plot, stat = stat_res))
} 

#' Calculate P values from a survival analysis
#' 
#' @param x result from either a chi square or cox regression
#' 
#' @return p value
#' @rdname pval-calc
#' @export
CalculatePValues.cox <- function(x) {
  logtest <- -2 * (x$loglik[1] - x$loglik[2])
  if (is.null(x$df)) 
    df <- sum(!is.na(coef))
  else df <- round(sum(x$df), 2)
  
  pval <- format.pval(pchisq(logtest, 
                             df, lower.tail = FALSE))
  return(pval)
}
