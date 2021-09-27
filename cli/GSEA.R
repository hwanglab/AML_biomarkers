#!/usr/bin/env Rscript
source("renv/activate.R")
library(argparse)

# parse args
parser <- ArgumentParser("Do GSEA")
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
  "--cores",
  "-c",
  help = "number of cores to use, 0 = all availible cores",
  default = 0
)
parser$add_argument(
  "--verbose",
  "-v",
  help = "should messages be printed? One of: DEBUG, INFO, WARN, ERROR",
  default = "INFO"
)

argv <- parser$parse_args(c("-i", "flt3_cebpa_stem_bc", "-v", "DEBUG"))

# load more libraries
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(DOSE)
  library(enrichplot)
  library(patchwork)
  library(glue)
  library(org.Hs.eg.db)
  library(future)
  library(tidyverse)
  library(furrr)
  library(log4r)
  library(here)
})

logger <- logger(threshold = argv$verbose)
if (argv$cores == 0) argv$cores <- availableCores()[[1]]
info(logger, glue("The number of workers is set to {argv$cores}"))

if (argv$dir == "") {
  output_path <- paste0("outs/", argv$id)
  plots_path <- paste0("plots/", argv$id)
} else {
  output_path <- paste0(parser$run_dir, "/outs/", argv$id)
  plots_path <- paste0(parser$run_dir, "/plots/", argv$id)
}

if (!dir.exists(here(plots_path))) {
  debug(logger, "Plots directory is being created")
  dir.create(here(plots_path))
}

if (!dir.exists(here(output_path))) {
  fatal(logger, "Output directory does not exist")
  quit()
}

markers <- read_tsv(here(output_path, "cluster_split_DE_master.tsv"), col_types = cols())

MakeGeneList <- function(df, entrez = FALSE) {
  gene_list <- df$avg_log2FC
  names(gene_list) <- as.character(df$gene)
  if (entrez) {
    names(gene_list) <- df$entrez
  }
  gene_list <- sort(gene_list, decreasing = TRUE)
  return(gene_list)
}

markers_fil <- filter(markers, is.finite(avg_log2FC)) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  mutate(entrez = mapIds(org.Hs.eg.db, gene, "ENTREZID", "SYMBOL"))

msigdb <- list(
  HALLMARK = "msigdb/gmt/h.all.v7.4.entrez.gmt",
  ONCO = "msigdb/gmt/c6.all.v7.4.entrez.gmt",
  REGULATORY = "msigdb/gmt/c3.tft.v7.4.entrez.gmt"
)
msigdb <- map(msigdb, read.gmt)

markers_fil_nest <- nest(markers_fil)

plan("multicore", workers = argv$cores)
furrr_options <- furrr_options(seed = 10000000, stdout = FALSE)
info(logger, "Preparing Gene List for GSEA")
markers_fil_nest <- mutate(
  markers_fil_nest,
  gene_list = map(data, MakeGeneList),
  gene_list_entrez = map(data, MakeGeneList, entrez = TRUE)
)

gene_list <- markers_fil_nest$gene_list
names(gene_list) <- markers_fil_nest$cluster

gene_list_entrez <- markers_fil_nest$gene_list_entrez
names(gene_list_entrez) <- markers_fil_nest$cluster

info(logger, "Running GSEA")
res <- list()
res[["GO"]] <- future_map(
  gene_list,
  ~ gseGO(
    geneList = .x,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    keyType = "SYMBOL",
    eps = 0,
    verbose = FALSE
  ),
  .options = furrr_options
)
res[["MKEGG"]] <- future_map(
  gene_list_entrez,
  ~ gseMKEGG(
    geneList = .x,
    organism = "hsa",
    eps = 0,
    verbose = FALSE
  ),
  .options = furrr_options
)
res[["KEGG"]] <- future_map(
  gene_list_entrez,
  ~ gseKEGG(
    geneList = .x,
    organism = "hsa",
    eps = 0,
    verbose = FALSE
  ),
  .options = furrr_options
)
res[["DO"]] <- future_map(
  gene_list_entrez,
  ~ gseDO(
    geneList = .x,
    eps = 0,
    verbose = FALSE
  ),
  .options = furrr_options
)
msigdb_res <- list()

for (set in names(msigdb)) {
  msigdb_res[[set]] <- future_map(
    gene_list_entrez,
    ~ GSEA(
      geneList = .x,
      TERM2GENE = msigdb[[set]],
      eps = 0,
      verbose = FALSE
    ),
    .options = furrr_options
  )
}

res <- append(res, msigdb_res)

plot_res <- transpose(res)
plot_res <- map(plot_res, ~ discard(.x, ~ !nrow(.x@result)))
dir.create(glue("{plots_path}/GSEA"))

for (cluster in names(plot_res)) {
  plots <- list()
  for (set in names(plot_res[[cluster]])) {
    p1 <- dotplot(plot_res[[cluster]][[set]], showCategory = 30)
    p2 <- cnetplot(
      plot_res[[cluster]][[set]],
      foldChange = gene_list_entrez[[cluster]]
    )
    edo <- pairwise_termsim(plot_res[[cluster]][[set]])
    p3 <- emapplot(edo)
    p4 <- upsetplot(edo)
    patch1 <- p1 / p4 + plot_layout(heights = c(2, 1))
    plots[[set]] <- p3 + patch1 +
      plot_annotation(
        title = glue("GSEA Plots for {cluster} | Pathway: {set}"),
        theme = theme(plot.title = element_text(size = 26))
      ) +
      plot_layout(widths = c(5, 2))
  }
  pdf(glue("{plots_path}/GSEA/{cluster}.pdf"), width = 35, height = 22)
  print(plots)
  graphics.off()
}

