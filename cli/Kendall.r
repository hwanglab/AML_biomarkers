#!/usr/bin/env Rscript
library(argparse)

# parse args
parser <- ArgumentParser("Assess Single Cell Gene Correlation")
parser$add_argument(
  "--dir",
  "-d",
  help = "path to run directory",
  default = ""
)
parser$add_argument(
  "--id",
  "-i",
  help = "ID to use for outputs, will read inputs from here"
)
parser$add_argument(
  "--verbose",
  "-v",
  help = "should messages be printed? One of: DEBUG, INFO, WARN, ERROR",
  default = "INFO"
)
parser$add_argument(
  "--genes",
  "-g",
  help = "genes to assess correlation of",
  nargs = "+"
)
argv <- parser$parse_args()

suppressPackageStartupMessages({
  library(Seurat)
  library(here)
  library(log4r)
  library(tidyverse)
  library(broom)
  library(glue)
})

logger <- logger(threshold = argv$verbose)

if (argv$dir == "") {
  output_path <- paste0("outs/", argv$id)
} else {
  output_path <- paste0(parser$run_dir, "/outs/", argv$id)
}

dir.create(glue("{output_path}/cor"))
# dir.create(glue("{plots_path}/cor"))
data_filename <- list.files(
  path = here(output_path, "cache/"),
  full.names = TRUE,
  pattern = "^seurat_dimred_"
)

if (length(data_filename) > 1) {
  fatal(logger, "There is more than one cached object")
  quit(status = 1)
}

seurat <- readRDS(data_filename)

DoKendall <- function(feat1, feat2, object) {
  data <- Seurat::FetchData(object, vars = c(feat1, feat2), slot = "data")
  cor_res <- cor.test(data[, 1], data[, 2], method = "kendall")
  return(cor_res)
}

info(logger, "Getting pairwise gene lists")
genes <- combn(argv$genes, 2) %>% t()
genes1 <- as.list(genes[, 1])
genes2 <- as.list(genes[, 2])

info(logger, "Doing correlation analysis")
cor_res <- map2(genes1, genes2, DoKendall, object = seurat)
names(cor_res) <- map2(genes1, genes2, paste, sep = "_")

info(logger, "Tidying statsical results")
res_tidy <- cor_res %>%
  map_df(broom::tidy) %>%
  mutate(var = names(cor_res)) %>%
  separate(var, into = c("gene1", "gene2")) %>%
  select(-method, -alternative, -statistic) %>%
  group_by(gene1) %>%
  arrange(p.value, .by_group = TRUE)

filename_feats <- glue_collapse(argv$genes, sep = "_", width = 10)
filename <- glue("{output_path}/cor/kendall_cor_{filename_feats}.tsv")

write_tsv(res_tidy, filename)
