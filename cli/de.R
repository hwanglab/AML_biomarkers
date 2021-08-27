#!/usr/bin/env Rscript
source("renv/activate.R")
library(argparse)

# parse args
parser <- ArgumentParser("Do Differential Expression")
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
  "--clusters",
  "-C",
  help = "clusters to use for DE",
  nargs = "+"
)
parser$add_argument(
  "--sig-level",
  "-s",
  help = "p-value threshold to use",
  default = 0.05,
  metavar = "P"
)
group <- parser$add_mutually_exclusive_group()
group$add_argument(
  "--B-mode",
  "-X",
  help = "Should B mode batch correction be applied?",
  action = "store_true"
)
group$add_argument(
  "--S-mode",
  "-S",
  help = "Should S mode batch correction be applied?",
  action = "store_true"
)

argv <- parser$parse_args()

# load more libraries
suppressPackageStartupMessages({
  library(here)
  library(future)
  library(Seurat)
  suppressWarnings(library(rsinglecell))
  library(pheatmap)
  library(readxl)
  library(tidyverse)
  library(furrr)
  library(log4r)
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

if (!dir.exists(here(plots_path))) {
  debug(logger, "Plots directory is being created")
  dir.create(here(plots_path))
}

if (!dir.exists(here(output_path))) {
  fatal(logger, "Output directory does not exist")
}

# read in data
data_filename <- list.files(
  path = here(output_path, "cache/"),
  full.names = TRUE,
  pattern = "^seurat_"
)

if (length(data_filename) > 1) fatal("There is more than one cached object")

debug(logger, paste0("Importing Data from: ", data_filename))
seurat <- readRDS(data_filename)

# set plan
if (argv$cores == 1) {
  plan("sequential")
} else {
  plan(tweak("multicore", workers = as.integer(argv$cores)))
}

info(logger, "Printing Plots")
WilcoxTestPossibly <- purrr::possibly(
  wilcox.test,
  otherwise = list(p.value = NA)
)

FindClusterFreq <- function(object, metadata, cluster, sort_by = NULL, .debug = FALSE) {
  # quote expressions
  md <- rlang::enquos(metadata)
  c <- rlang::enquo(cluster)
  if (.debug == TRUE) browser()
  # select meta.data columns from seurat object and calculate freq
  if (is.null(sort_by)) {
    freq <- object %>%
      dplyr::group_by(dplyr::across(.cols = all_of(c(metadata, cluster)))) %>%
      summarise(n = n()) %>%
      dplyr::group_by(dplyr::across(.cols = all_of(c(metadata)))) %>%
      mutate(freq = n / sum(n) * 100)
  } else {
    sb <- rlang::enquo(sort_by)
    freq <- object %>%
      dplyr::group_by(dplyr::across(.cols = all_of(c(metadata, cluster)))) %>%
      summarise(n = n(), .groups = "keep") %>%
      mutate(freq = n / sum(n) * 100) %>%
      dplyr::arrange(!!sb)
  }
  return(freq)
}
freq <- FindClusterFreq(
  seurat[[]],
  c("patient_id", "prognosis"),
  "clusters"
)

freq_stat <- freq %>%
  group_by(clusters) %>%
  summarise(
    pval = WilcoxTestPossibly(
      x = freq,
      y = as.numeric(as.factor(prognosis))
    )$p.value,
    .groups = "keep"
  )

wilcox_clusters <- freq_stat %>%
  filter(pval <= argv$sig_level) %>%
  pull(clusters)

print_clusters <- glue_collapse(wilcox_clusters, sep = ", ", last = " and ")
debug(logger, glue("The following clusters are selected: {print_clusters}"))

bad_cells <- WhichCells(seurat, idents = wilcox_clusters, invert = TRUE)
Idents(seurat, cells = bad_cells) <- "Not Significant"

info(logger, "Printing Dimension Reductions")
pdf(here(plots_path, "UMAP.pdf"), onefile = TRUE, width = 10)
DimPlot(seurat, reduction = "umap", group.by = "clusters")
DimPlot(seurat, reduction = "umap", split.by = "prognosis")
graphics.off()

info(logger, "Printing Cluster Bar Plots")
bar_plot <- ggplot(
  data = freq,
  mapping = aes(x = patient_id, y = freq, fill = prognosis)
) +
  geom_col() +
  theme_classic() +
  scale_fill_viridis_d() +
  facet_wrap(~clusters) +
  scale_x_discrete(labels = NULL, breaks = NULL)

bar_plot2 <- freq %>%
  filter(clusters %in% wilcox_clusters) %>%
  ggplot(mapping = aes(x = patient_id, y = freq, fill = prognosis)) +
  geom_col() +
  theme_classic() +
  scale_fill_viridis_d() +
  facet_wrap(~clusters) +
  scale_x_discrete(labels = NULL, breaks = NULL)

freq2 <- freq %>%
  group_by(clusters, prognosis) %>%
  summarise(x_bar = mean(freq, na.rm = TRUE), sd = sd(freq))

bar_plot3 <- freq2 %>%
  filter(clusters %in% wilcox_clusters) %>%
  ggplot(mapping = aes(x = prognosis, y = x_bar, fill = prognosis)) +
  geom_col() +
  theme_classic() +
  scale_fill_viridis_d() +
  facet_wrap(~clusters) +
  scale_x_discrete(labels = NULL, breaks = NULL)

bar_plot4 <- freq2 %>%
  filter(clusters %in% wilcox_clusters) %>%
  ggplot(mapping = aes(x = prognosis, y = x_bar, fill = prognosis)) +
  geom_col() +
  theme_classic() +
  scale_fill_viridis_d() +
  facet_wrap(~clusters) +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  geom_errorbar(
    mapping = aes(ymin = x_bar - sd, ymax = x_bar + sd),
    width = 0.2,
    position = position_dodge(0.9)
  )

pdf(here(plots_path, "cluster_bar_plot.pdf"))
print(bar_plot)
print(bar_plot2)
print(bar_plot3)
print(bar_plot4)
graphics.off()

# info(logger, "Printing single-cell Heatmap")
# pdf(here(plots_path, "feature_heatmap.pdf"), width = 18)
# print(DoHeatmap(seurat, group.by = "clusters"))
# graphics.off()

if (!file.exists(here(output_path, "cluster_differential_expression.tsv"))) {
  info(logger, "Doing Differential Expression")
  markers <- FindAllMarkers(seurat, test.use = "MAST")
  write_tsv(
    markers,
    file = here(output_path, "cluster_differential_expression.tsv")
  )
} else {
  info(logger, "Differential Expression already done")
  markers <- read_tsv(
    file = here(output_path, "cluster_differential_expression.tsv"),
    col_types = cols()
  )
}

ReturnGEP <- function(dataset) {
  fname <- here(
    output_path,
    cibersort_results_dir,
    glue("CIBERSORTxGEP_{dataset}_SM_GEPs_Filtered.txt")
  )
  gep <- readr::read_tsv(fname, col_types = readr::cols())
  return(gep)
}

datasets <- c("target_data.txt", "tcga_data.txt", "beat_aml.txt")

if (argv$B_mode) {
  cibersort_results_dir <- "cibersort_results_B"
} else {
  if (argv$S_mode) {
    cibersort_results_dir <- "cibersort_results_S"
  } else {
    cibersort_results_dir <- "cibersort_results"
  }
}

gene_names_for_de <- list(0)
for (set in datasets) {
  do_heatmap <- TRUE
  # Probably not the most efficent
  tryCatch(ReturnGEP(set),
    error = function(e) {
      warn(logger, glue("CIBERSORTxGEP output for {set} not found"))
      do_heatmap <<- FALSE
    }
  )
  if (do_heatmap) {
    gep <- ReturnGEP(set)
    gep <- column_to_rownames(gep, var = "GeneSymbol")
    gep[is.na(gep)] <- 0
    gep_sd <- apply(gep, 1, sd)
    q3 <- quantile(gep_sd)[[4]]
    gep_bar <- mean(gep_sd)

    gep_filter <- gep[!(gep_sd <= gep_bar), ]
    gene_names_for_de[[set]] <- rownames(gep_filter)
  }
}
gene_names_for_de <- unique(unlist(gene_names_for_de))

if (!is_empty(gene_names_for_de)) {
  info(logger, "Scaling Genes from CIBERSORTxGEP")
  debug(logger, glue("Getting Pearson Residuals for {length(gene_names_for_de)} genes"))
  seurat <- GetResidual(seurat, gene_names_for_de, verbose = FALSE)
  
  debug(logger, "Creating Heatmap")
  pdf(
    here(plots_path, "feature_heatmap_CIBERSORTxGEP_filtered.pdf"),
    width = 18
  )
  print(DoHeatmap(seurat, group.by = "clusters", features = gene_names_for_de))
  graphics.off()
}

genes <- markers %>%
  group_by(cluster) %>%
  slice_min(p_val_adj, n = 50) %>%
  pull(gene)

if (!is_empty(genes)) {
  info(logger, "Scaling Genes from DE")
  seurat <- GetResidual(seurat, genes, verbose = FALSE)
  info(logger, "Printing DE filtered single-cell Heatmap")
  pdf(here(plots_path, "feature_heatmap_DE_genes.pdf"), width = 18)
  print(DoHeatmap(seurat, group.by = "clusters", features = genes))
  graphics.off()
}

if (length(argv$clusters) != 0) {
  if (length(argv$clusters) == 1) {
    error(logger, "You must have at least 2 clusters selected")
    quit()
  }
  cells <- WhichCells(seurat, idents = argv$clusters)
  seurat_selected <- seurat[, cells]
  markers <- FindAllMarkers(seurat_selected, test.use = "MAST")
  clusters <- glue_collapse(clusters, sep = ", ", final = " and ")
  fname <- glue("cluster_differential_expression_({clusters}).tsv")
  write_tsv(markers, file = here(output_path, fname))
}

source(here("lib/WriteInvocation.R"))
WriteInvocation(argv, output_path = here(output_path, "invocation"))
