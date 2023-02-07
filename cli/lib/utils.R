
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
      paste(pkgs[!package.installed], collapse = ", "),
      ". Please install"
    )
  }
  invisible(x = package.installed)
}

#' Not in operator
#'
#' @param lhs stuff
#' @param rhs stuff
#'

`%!in%` <- function(lhs, rhs) {
  match(lhs, rhs, nomatch = 0) == 0
}

#' Parse information for survival
#' @param s a string see docs for formatting
ReturnNamesFromString <- function(s) {
  data_info <- c(
    "TARGET:efs=Event Free Survival Time in Days:wbc=wbc_at_diagnosis:status=First Event:not_event=censored", # nolint
    "TCGA:efs=days_to_death:status=vital_status:event=Dead",
    "BeatAML:efs=OS_DAYS:status=status:event=1"
  )

  if (length(s) != 0) data_info <- c(data_info, s)
  split <- stringr::str_split(data_info, ":")
  names <- purrr::map_chr(split, ~ .x[[1]])
  split1 <- purrr::map(split, ~ .x[-1])

  split2 <- purrr::map(split1, stringr::str_split, "=")
  names(split2) <- names
  purrr::map(
    split2,
    ~ purrr::map(.x, ~ rlang::set_names(.x[[2]], .x[[1]])) %>% unlist()
  ) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(set = names)
}

ReturnNamesFromJSON <- function(s, training) {
  if (!file.exists(s)) {
    error(logger, glue("File does not exist: {s}"))
    quit(status = 1)
  }
  p <- jsonlite::fromJSON(s)
  if (!training) p <- map_dfr(p, as_tibble, .id = "set")
  return(p)
}

ParseNames <- function(s, training = FALSE) {
  if (length(s) == 0) {
    return(ReturnNamesFromString(s))
  } else {
    if (tools::file_ext(s) == "json") {
      return(ReturnNamesFromJSON(s, training))
    } else {
      return(ReturnNamesFromString(s))
    }
  }
}


#' Get the filename of a saved model
#' @param type a string used to identify the model
GetSavedModelFilename <- function(type) {
  if (type == "nn_regression") {
    here(output_path, "saved_models", glue("{type}_model.h5"))
  } else {
    here(output_path, "saved_models", glue("{type}_model.rds"))
  }
}

#' Check if model has been saved
#' @inheritParams GetSavedModelFilename
CheckSavedModel <- function(type) {
  model_save <- GetSavedModelFilename(type)
  file.exists(model_save) || dir.exists(model_save)
}

#' Load a saved model from disk
#' @inheritParams GetSavedModelFilename
LoadSavedModel <- function(type) {
  model_save <- GetSavedModelFilename(type)
  if (type == "nn_regression") {
    backend_k <- keras::backend()
    backend_k$clear_session()
    model <- keras::load_model_hdf5(model_save)
  } else {
    model <- readRDS(model_save)
  }
  return(model)
}

#' Load Annotated Deconvolution Matrix
#' @param output_path output directory where information is stored
#' @param logger a logger object for printing progress
LoadAnnotatedMatrix <- function(output_path, logger) {
  data_filename <- here(output_path, glue("cache/clinical_deconvoluted.rds"))
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
  return(data)
}

#' Stop Execution if Output directory does not exist
StopIfOutputDirNotExist <- function(output_path) {
  if (!dir.exists(here(output_path))) {
    fatal(logger, "Output directory does not exist")
    quit(status = 1)
  }
}

#' Get output directory from argument parser
#' @param argv and ArgumentParser object or list with indexes for `dir` and `id`
PrepareOutDir <- function(argv) {
  if (argv$dir == "") {
    output_path <- paste0("outs/", argv$id)
  } else {
    output_path <- paste0(parser$dir, "/outs/", argv$id)
  }
  return(output_path)
}

#' Get and validate parameters for training data
GetAndValidateTrainingParams <- function(output_path, argv, data) {
  params <- ParseNames(argv$training, training = TRUE)
  file <- here(glue("{output_path}/saved_models/training_hash.rds"))

  train_vars <- params$use_for_training
  print(train_vars)
  params$use_for_training <- NULL
  params <- list(
    df = data,
    cols = params,
    vars = train_vars,
    time_unit = "months"
  )
  hash <- digest::digest(params)
  if (!file.exists(file)) {
    saveRDS(hash, file)
    return(params)
  }
  old_hash <- readRDS(file)

  if (old_hash != hash) {
    error(logger, "You have supplied a different set of paramters than were supplied previously.")
    info(logger, "You should either run with the old params or select a new test id")
    quit(status = 1)
  }
  return(params)
}