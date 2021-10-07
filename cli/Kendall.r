#!/usr/bin/env Rscript
source("renv/activate.R")
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
  "--batch-correct",
  "-X",
  help = "Should intetegraed data be used from Seurat",
  action = "store_true"
)
parser$add_argument(
  "--genes",
  "-g",
  help = "genes to assess correlation of",
  nargs = "+"
)

suppressPackageStartupMessages({
  library(Seurat)
  library(here)
  library(log4r)
  library(tidyverse)
  library(glue)
})



logger <- logger(threshold = argv$verbose)

if (argv$dir == "") {
  output_path <- paste0("outs/", argv$id)
  plots_path <- paste0("plots/", argv$id)
} else {
  output_path <- paste0(parser$run_dir, "/outs/", argv$id)
  plots_path <- paste0(parser$run_dir, "/plots/", argv$id)
}

create.dir(glue("{output_path}/cor"))
# create.dir(glue("{plots_path}/cor"))
data_filename <- list.files(
  path = here(output_path, "cache/"),
  full.names = TRUE,
  pattern = "^seurat_dimred_"
)

if (length(data_filename) > 1) {
  fatal(logger, "There is more than one cached object")
  quit()
}

seurat <- readRDS(data_filename)

p <- list()
res <- list()
for (j in seq_along(features)) {
  data_feature <- FetchData(seurat, vars = c(features[[j]], "BCL2"))

  res[[j]] <- cor.test(data_feature[, 1], data_feature[, 2], method = "kendall")

  data2 <- data_feature %>% as.matrix() %>% t() %>% as.data.frame()

  p[[[j]] <- ggplot(data = data2, mapping = aes_string(x = features[[j]], y = "BCL2")) +
    geom_point() +
    theme_classic()
  }
names(res) <- features
names(p) <- features

DoKendall <- function(feat1, feat2, object) {
    data <- Seurat::FetchData(object, vars = c(feat1, feat2), slot = "data")
    cor_res <- cor.test(data[, 1], data[, 2], method = "kendall")
    return(cor_res)
}

y <- combn(x, 2) %>% t()
feats1 <- as.list(y[, 1])
feats2 <- as.list(y[, 2])

cor_res <- map2(feats1, feats2, DoKendall, object = seurat)
names(cor_res) <- map2(feats1, feats2, paste, sep = "_")

res_tidy <- cor_res %>%
  map_df(broom::tidy) %>%
  mutate(var = names(cor_res)) %>%
  separate(var, into = c("gene1", "gene2")) %>%
  select(-method, -alternative, -statistic) %>%
  group_by(gene1) %>%
  arrange(p.value, .by_group = TRUE)

filename_feats <- glue_collapse(features, sep = "_", width = 10)
filename <- glue("{output_path}/cor/kendall_cor_{filename_feats}.tsv")

write_tsv(res_tidy, filename)