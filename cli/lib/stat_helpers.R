#' Prepare a cox-proportonal hazards model and forest plot
#' @param d data for modeling
#' @param formula a formula for coxph
DoCox <- function(d, formula) {
    fit <- survival::coxph(formula, data = d)
    title <- glue::glue("Liklihood Ratio = {AngrilyExtractLiklihoodRatio(fit)}")
    plot <- ggforest(fit, data = d, main = title)
    return(list(fit = fit, plot = plot))
}

#' Get Liklihood ratio from a coxph model
#' @param fit a coxph fit object
# angry because: https://twitter.com/RobertSchauner/status/1516419829053464583
AngrilyExtractLiklihoodRatio <- function(fit) {
    if (!("coxph" %in% class(fit))) {
        return(invisible(NULL))
    }
    return(-2 * (fit$loglik[1] - fit$loglik[2]))
}
