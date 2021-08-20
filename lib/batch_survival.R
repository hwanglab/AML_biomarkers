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
            debug(logger, paste0("Using ", avail_cols[[st]], " on ", names(data)[[dataset]]))
            x[["status"]] <- if_else(x[[avail_cols[[st]]]] == avail_events[[st]], 1, 0)
            x[["score"]] <- UseLASSOModelCoefs(x, coef(lasso_model))
            x[["score_bin"]] <- if_else(x$score >= median(x$score, na.rm = TRUE), "High", "Low")
            res <- list(0)

            for (time_axis in seq_along(times)) {
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
            names(res) <- times
            res <- purrr::discard(res, ~ all(is.na(.x)))

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
        plots_final[[dataset]] <- plots
        stats_final[[dataset]] <- reduce(stats, bind_rows)
        score_final[[dataset]] <- x[["score"]]
    }
    res1 <- stats_final %>% keep(is.data.frame) %>% reduce(bind_rows)
    res2 <- plots_final %>% unlist(recursive = FALSE) %>% unlist(recursive = FALSE)
    return(list(stats = res1, plots = res2, scores = score_final))
}
