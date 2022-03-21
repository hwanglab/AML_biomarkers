
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
#'
ReturnNamesFromString <- function() {
    if (length(argv$info) != 0) data_info <- c(data_info, argv$info)
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

ReturnNamesFromJSON <- function() {
    if (!file.exists(argv$info)) {
        error(logger, glue("File does not exist: {argv$info}"))
        quit(status = 1)
    }
    jsonlite::fromJSON(argv$info) %>% map_dfr(as_tibble, .id = "set")
}

ParseNames <- function() {
    if (length(argv$info) == 0) {
        return(ReturnNamesFromString())
    } else {
        if (tools::file_ext(argv$info) == "json") {
            return(ReturnNamesFromJSON())
        } else {
            return(ReturnNamesFromString())
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
