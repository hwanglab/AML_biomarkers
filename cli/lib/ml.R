
#' Clean data to prepare for training and prediction
#'
#' @param df data to clean
#' @param wbc column with % WBCs as string
#' @param efs column with survival as string
#' @param status column with status as string
#' @param time_unit string with unit of time to use
CleanData <- function(df, wbc = NULL, efs = NULL, status, time_unit = "days") {
  clean <- janitor::clean_names(df)

  wbc <- janitor::make_clean_names(wbc)
  efs <- janitor::make_clean_names(efs)
  status <- janitor::make_clean_names(status)

  namekey <- c("original_time", "status")
  names(namekey) <- c(efs, status)

  clean <- plyr::rename(clean, namekey, warn_missing = TRUE) %>%
    mutate(
      days = lubridate::days(as.integer(original_time)),
      time = lubridate::time_length(days, unit = time_unit) %>% as.numeric(),
      days = original_time
    ) %>%
    select(starts_with("cluster"), time, status, days)

  # Suppress "Removed n rows containing missing values." warnings
  suppressWarnings(
    res <- clean %>%
      as_tibble(.name_repair = "minimal") %>%
      setNames(colnames(.)) %>%
      remove_missing()
  )
  return(res)
}

#' Subset data with good features
#' @param df a datafame
PrepareDataForML <- function(df, categorical = FALSE) {
  data <- dplyr::select(df, starts_with("cluster"))

  if (categorical) {
    labs <- pull(df, status)
  } else {
    labs <- pull(df, days)
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
  res <- dplyr::select(data, status, time)
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

SetupTARGETData <- function(data, surv_data) {
  res <- data %>%
    mutate(
      AR = if_else(`FLT3/ITD allelic ratio` > 0.4, "AR+", "AR-"),
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
    ) %>%
    select(any_of(mutation_cols), AR) %>%
    bind_cols(surv_data)
  return(res)
}
