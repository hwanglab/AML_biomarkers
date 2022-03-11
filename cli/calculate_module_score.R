#!/usr/bin/env Rscript
library(argparse)

# parse args
parser <- ArgumentParser("Calculate Module Scores across Multiple IDs")
parser$add_argument(
  "--dir",
  "-d",
  help = "path to run directory",
  default = ""
)
parser$add_argument(
  "--ids",
  "--id",
  "-i",
  help = "ID to read inputs, output will be under outs in dir",
  nargs = "+"
)
parser$add_argument(
  "--verbose",
  "-v",
  help = "should messages be printed?",
  default = "INFO",
  choices = c("DEBUG", "INFO", "WARN", "ERROR")
)
parser$add_argument(
  "--cores",
  "-c",
  help = "how many cores should be used",
  default = 1
)
parser$add_argument(
  "--slurm",
  "-S",
  help = "should batchtools be used with slurm",
  action = "store_true"
)
parser$add_argument(
  "--genes",
  "-G",
  help = "path to gene list file as gtf"
)


argv <- parser$parse_args(c("-i", "flt3_cebpa_stem_bc", "flt3_cebpa_nonstem_bc", "-G", "lsc_genesets.gmt"))

# load more libraries
suppressPackageStartupMessages({
  library(here)
  library(future)
  library(Seurat)
  library(readxl)
  library(qusage)
  library(tidyverse)
  library(furrr)
  library(log4r)
  library(glue)
  library(xfun)
})

logger <- logger(threshold = argv$verbose)


if (argv$dir == "") {
  output_path <- paste0("outs")
  plots_path <- paste0("plots")
} else {
  output_path <- paste0(parser$run_dir, "/outs")
  plots_path <- paste0(parser$run_dir, "/plots")
}

if (!dir.exists(here(plots_path))) {
  debug(logger, "Plots directory is being created")
  dir.create(here(plots_path))
}

if (!dir.exists(here(output_path))) {
  fatal(logger, "Output directory does not exist")
}

# read in data

inputs <- map(argv$ids, ~ glue("{output_path}/{.x}/cache"))
data_filename <- map(
  inputs,
  ~ list.files(
    path = here(.x),
    full.names = TRUE,
    pattern = "^seurat_dim"
  )
)

seurats <- lapply(data_filename, readRDS)
seurat <- reduce(seurats, merge)


seurat <- cache_rds(
  SCTransform(seurat, conserve.memory = TRUE),
  file = "seurat_merged_normalized.rds",
  hash = list(seurat[["RNA"]])
)

# Suppress readLines warning about incomplete final line --its ok here
sets <- suppressWarnings(read.gmt(argv$genes))

seurat <- AddModuleScore(seurat, features = sets, search = TRUE, name = "ModuleScore")

md <- seurat[[]] %>%
  select(starts_with("ModuleScore")) %>%
  set_names(names(sets)) %>%
  mutate(group = pull(seurat[["stemness"]], 1)) %>%
  as_tibble()

library(patchwork)

fatal(logger, "Script unfinished")
quit(status = 1)
for (i in seq_along(names(set))) {
  plot[[i]] <- ggplot(data = md, mapping = aes_string(y = names(sets)[[i]], x = "group")) +
    geom_violin()
}

ggplot(data = md, mapping = aes_string(y = names(sets)[[1]], x = "group")) +
  geom_violin() +
  theme_classic()
