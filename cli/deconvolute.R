#!/usr/bin/env Rscript
library(argparse)

# parse args
parser <- ArgumentParser("Deconvolute Samples")
parser$add_argument(
  "--dir",
  "-d",
  help = "path to run directory",
  default = ""
)
parser$add_argument(
  "--ref",
  "-r",
  help = "path to CIBERSORTx GEP",
  default = ""
)
parser$add_argument(
  "--num-cells",
  "-n",
  help = "number of cells to use for ref creation",
  default = 100
)
parser$add_argument(
  "--id",
  "-i",
  help = "ID to use for outputs"
)
parser$add_argument(
  "--cores",
  "-c",
  help = "number of cores to use",
  default = 1
)
parser$add_argument(
  "--verbose",
  "-v",
  help = "should messages be printed? One of: DEBUG, INFO, WARN, ERROR",
  default = "INFO"
)
parser$add_argument(
  "--plot-similarity",
  "-s",
  help = "should the reference simelarity be plotted?",
  action = "store_true"
)

argv <- parser$parse_args()

# load more libraries
suppressWarnings({
  suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratDisk)
    library(here)
    library(xfun)
    library(granulator)
    library(log4r)
    library(rsinglecell)
    library(tidyverse)
  })
})

if (argv$dir == "") {
  output_path <- paste0("outs/", argv$id)
  plots_path <- paste0("plots/", argv$id)
} else {
  output_path <- paste0(parser$run_dir, "/outs/", argv$id)
  plots_path <- paste0(parser$run_dir, "/plots/", argv$id)
}
if (!dir.exists(here(output_path))) {
  fatal(logger, "Output directory does not exist")
}

logger <- logger(threshold = argv$verbose)

data_filename <- list.files(
  path = here(output_path, "cache"),
  full.names = TRUE,
  pattern = "^seurat_"
)

if (length(data_filename) > 1) fatal("There is more than one cached object")

debug(logger, paste0("Importing Data from: ", data_filename))
diagnosis <- readRDS(data_filename)

info(logger, "Making References")

refs <- cache_rds(
  {
    debug(logger, "Starting Expression Testing")
    FindAllMarkers <- quietly(FindAllMarkers)

    ref <- FindAllMarkers(
      diagnosis,
      max.cells.per.ident = argv$num_cells,
      test.use = "MAST"
    )$result
    info(logger, "Expression Testing Done")
    ref2 <- ref[ref$p_val_adj <= 0.05 & abs(ref$avg_log2FC) >= 0.25, "gene"]
    debug(logger, "Reference Filtered")
    ref3 <- AverageExpression(diagnosis,
      assays = "SCT",
      features = ref2,
      group.by = "clusters"
    )$SCT

    debug(logger, "Expression Averaged")
    colnames(ref3) <- colnames(ref3) %>% paste0("cluster", .)

    gep_filename <- list.files(
      path = here(output_path),
      full.names = TRUE,
      pattern = "^CIBERSORTx_Job"
    )

    if (is.null(gep_filename)) {
      error(logger, "Cannot find CIBERSORTx GEP")
    }

    debug(logger, paste0("Loading CIBERSORTx GEP from:", gep_filename))
    ciber_gep <- read_tsv(here(gep_filename), col_types = cols()) %>%
      column_to_rownames(var = "GeneSymbol") %>%
      as.matrix()

    list(CIBERSORTx = ciber_gep, custom = ref3)
  },
  file = paste0("decon_refs.rds"),
  dir = paste0(output_path, "/cache/"),
  hash = list(diagnosis[["clusters"]])
)

if (argv$plot_similarity) {
  info(logger, "Plotting similarity")
  sim_plot <- plot_similarity(refs)
} else {
  info(logger, "Skipping similarity plotting")
}

info(logger, "Reading in Bulk Data")
target_data <- read_tsv(here("cibersort_in/target_data.txt"), col_types = cols()) %>%
  column_to_rownames(var = "GeneSymbol") %>%
  as.matrix()

tcga_data <- read_tsv(here("cibersort_in/tcga_data.txt"), col_types = cols()) %>%
  column_to_rownames(var = "GeneSymbol") %>%
  as.matrix()

beatAML_data <- read_tsv(here("cibersort_in/beat_aml.txt"), col_types = cols()) %>%
  column_to_rownames(var = "GeneSymbol") %>%
  as.matrix()

methods <- get_decon_methods()[get_decon_methods() %!in% c("svr")] %>%
  as.vector()

datasets <- list(TARGET = target_data, TCGA = tcga_data, BeatAML = beatAML_data)
datasets <- map(
  datasets,
  ~ .x %>%
    as.data.frame() %>%
    distinct() %>%
    as.matrix()
)

info(logger, "Deconvoluting Samples")
deconvoluted_samples <- cache_rds(
  map(
    datasets,
    ~ suppressMessages(
      deconvolute(
        .x,
        sigMatrix = refs,
        methods = methods,
        use_cores = argv$cores
      )
    )
  ),
  hash = list(datasets, refs),
  file = "deconvoluted_samples.rds",
  dir = paste0(output_path, "/cache/")
)

### Make + Print Plots ----

info(logger, "Preparing Plots")
decon_plots <- map(deconvoluted_samples,
  plot_deconvolute,
  scale = TRUE,
  labels = FALSE,
  markers = FALSE
)
cor_res <- map(deconvoluted_samples, correlate)
cor_plots <- map(cor_res,
  plot_correlate,
  method = "heatmap",
  legend = TRUE
)

decon_plots2 <- map(decon_plots, .f = function(ggplot) {
  data <- ggplot$data
  data <- filter(data, celltype == "cluster0")
  ggplot$data <- data
  return(ggplot + facet_wrap(~model))
})

pdf(file = here(plots_path, "decon_plots.pdf"), width = 12, height = 18)
sim_plot
decon_plots
decon_plots2
cor_plots
graphics.off()

source(here("lib/WriteInvocation.R"))
WriteInvocation(argv, output_path = here(output_path, "invocation"))
