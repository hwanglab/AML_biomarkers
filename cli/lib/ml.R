
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
PrepareDataForML <- function(df, categorical = FALSE) {
  data <- dplyr::select(df, starts_with("cluster"))

  if (categorical) {
    labs <- pull(df, status)
  } else {
    labs <- pull(df, efs_days)
  }

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