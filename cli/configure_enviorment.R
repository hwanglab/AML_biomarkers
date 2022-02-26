#!/usr/bin/env Rscript

DownloadPackage <- function(package) {
    install.packages(
        package,
        repos = "https://cloud.r-project.org",
        quiet = TRUE
    )
}

suppressPackageStartupMessages({
    if (!require("argparse", quietly = TRUE)) DownloadPackage("argparse")
    if (!require("log4r", quietly = TRUE)) DownloadPackage("log4r")
    if (!require("renv", quietly = TRUE)) DownloadPackage("renv")
    library(argparse)
    library(log4r)
})

parser <- ArgumentParser("Setup R Package Enviorment")

argv <- parser$parse_args()

logger <- logger(threshold = "INFO")

info(logger, "Checking if renv has already been used")

tryCatch(
    consent <- renv:::renv_consent_check(),
    error = function(e) {
        warn(logger, "Providing consent to renv")
        renv::consent(provided = TRUE)
    }
)

info(logger, "Now installing packages. This could take some time")
renv::restore(prompt = FALSE)

info(logger, "Packages installed successfully!")
info(logger, "Activating project")
renv::activate()

info(logger, "R enviorment set up!")
