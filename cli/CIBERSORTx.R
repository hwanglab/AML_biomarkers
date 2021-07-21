#!/usr/bin/env -S Rscript --no-save --quiet
source("renv/activate.R")
library(argparser)

# parse args
parser <- arg_parser("Deconvolute Samples")
parser <- add_argument(
  parser, "--dir",
  short = "-d",
  help = "path to run directory",
  default = ""
)
parser <- add_argument(
  parser, "--ref",
  short = "-r",
  help = "path to CIBERSORTx GEP",
  default = ""
)
parser <- add_argument(
  parser,
  "--id",
  short = "-i",
  help = "ID to use for outputs"
)
parser <- add_argument(
  parser,
  "--cores",
  short = "-c",
  help = "number of cores to use",
  default = 8
)
parser <- add_argument(
  parser,
  "--verbose",
  short = "-v",
  help = "should messages be printed? One of: DEBUG, INFO, WARN, ERROR",
  default = "INFO"
)
parser <- add_argument(
  parser,
  "--use-slurm",
  short = "-S",
  help = "should slurm be used for parallelization?",
  flag = TRUE
)
argv <- parse_args(parser)

suppressPackageStartupMessages({
  library(here)
  library(log4r)
})

logger <- logger(argv$verbose)

use_singularity <- FALSE
use_singularity <- tryCatch(
  {
    system("docker")
  },
  warning = function(w) {
    TRUE
  }
)

if (!is.logical(use_singularity)) use_singularity <- FALSE

debug(logger, paste0("Docker Tested with this final value for use_singularity: ", use_singularity))

if (argv$dir == "") {
  output_path <- paste0("outs/", argv$id)
} else {
  output_path <- paste0(parser$run_dir, "/outs/", argv$id)
}
if (!dir.exists(here(output_path))) {
  fatal(logger, "Output directory does not exist")
}

if (!file.exists(here(output_path, "cibersort_ref_input.txt"))) {
  info(logger, "CIBERSORTx Prepared files not detected.")

  sys_cmd <- paste(
    "./cli/make_CIBERSORT_ref.R",
    "-i", argv$id,
    "-d", argv$dir,
    "-v", argv$verbose
  )
  system(sys_cmd)
}

debug(logger, "Copying Data")
files <- list.files(here("cibersort_in"), full.names = TRUE)
file.copy(files, here(output_path))

info(logger, "Preparing commands for CIBERSORTx")


if (use_singularity) {
  info(logger, "Docker is not installed, falling back to singularity.")
  debug(logger, "Testing singularity")
  tryCatch(
    {
      system("singularity")
    },
    warning = function(w) {
      error(logger, "Singularity is not installed")
    }
  )
  cmd_start <- paste(
    "singularity run --pwd /src",
    "-B", paste0(here(output_path), ":/src/data"),
    "-B", paste0(here(output_path, "cibersort_results"), ":/src/outdir"),
    "docker://cibersortx/fractions"
  )
} else {
  info(logger, "Docker detected.")
  cmd_start <- paste(
    "docker run",
    "-v", paste0(here(output_path, "cibersort_in"), ":/src/data"),
    "-v", paste0(here(output_path, "cibersort_results"), ":/src/outdir"),
    "cibersortx/fractions"
  )
}

debug(logger, paste0("Using the username: ", Sys.getenv("EMAIL")))
debug(logger, paste0("Using the token: ", Sys.getenv("TOKEN")))

base_cmd <- paste(
  cmd_start,
  "--username", Sys.getenv("EMAIL"),
  "--verbose", "FALSE",
  "--token", Sys.getenv("TOKEN"),
  "--single_cell", "TRUE",
  "--outdir", here(output_path, "cibersort_results")
)

if (argv$verbose != "DEBUG") {
  debug_redir <- "1> /dev/null"
} else {
  debug_redir <- ""
}

sys_cmd <- paste(
  base_cmd,
  "--refsample", "cibersort_ref_input.txt",
  debug_redir
)
ref_filename <- "CIBERSORTx_cibersort_ref_input_inferred_refsample.txt"

if (!file.exists(here(output_path, "cibersort_results", ref_filename))) {
  info(logger, "Making CIBERSORTx Reference")
  system(sys_cmd)
} else {
  info(logger, "Reference Found. Skipping Reference Creation")
}

file.copy(here(output_path, "cibersort_results", ref_filename), here(output_path, "cibersort_in"))

info(logger, "Running CIBERSORTx")

sys_cmd <- function(mixture) {
  paste(
    base_cmd,
    "--sigmatrix", ref_filename,
    "--mixture", mixture,
    debug_redir
  )
}

if (argv$use_slurm) {
  if (!require("slurmR")) {
    slurm_install <- tryCatch(devtools::install_github("USCbiostats/slurmR"),
      error = function(e) {
        info(logger, "slurmR not installed. Not using slurm")
        NULL
      }
    )
    if (is.null(slurm_install)) {
      argv$use_slurm <- FALSE
    }
  }
}

if (argv$use_slurm) {
  lapply_cus <- function(X, FUN) {
    Slurm_lapply(X, FUN, plan = "collect", njobs = length(X))
  }
} else {
  lapply_cus <- function(X, FUN) lapply(X = X, FUN = FUN)
}

files_list <- list("beat_aml_1_sam_test.txt")

cibersort_results <- lapply_cus(X = files_list, FUN = function(x) {
  system(sys_cmd(x))
})

source(here("lib/WriteInvocation.R"))
WriteInvocation(argv, output_path = here(output_path, "invocation"))