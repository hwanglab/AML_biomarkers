#!/usr/bin/env Rscript
library(argparser)

# parse args
parser <- arg_parser("Do Survival Analysis")
parser <- add_argument(
  parser, "--dir",
  short = "-d",
  help = "path to run directory",
  default = ""
)
parser <- add_argument(
  parser,
  "--id",
  short = "-i",
  help = "ID to use for outputs",
  default = "test"
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
  "--cores",
  short = "-c",
  help = "number of cores",
  default = 12
)

argv <- parse_args(parser)

# load more libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(here)
  library(log4r)
  library(survival)
  library(survminer)
  library(scorecard)
  library(glmnet)
  library(rsinglecell)
  library(tidyverse)
})

source(here("lib/functions.R"))

logger <- logger(threshold = argv$verbose)

source(here("cli/lib/utils.R"))
output_path <- PrepareOutDir(argv)
StopIfOutputDirNotExist(output_path)

plots_path <- here(output_path, "plots")

# read in data
data_filename <- list.files(
  path = here(output_path, "cache/"),
  full.names = TRUE,
  pattern = "^seurat_"
)

if (length(data_filename) > 1) fatal("There is more than one cached object")

debug(logger, paste0("Importing Data from: ", data_filename))
diagnosis <- readRDS(data_filename)

tmp <- LoadH5Seurat(
  file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat"),
  assays = c("RNA")
)
DefaultAssay(tmp) <- "RNA"

debug(logger, "Subsetting Seurat Object")

tmp <- tmp[, Cells(diagnosis)]
tmp <- NormalizeData(tmp)

Idents(diagnosis) %>%
  as_tibble(rownames = "Cell") %>%
  rename(cell_type = value) %>%
  write_tsv(here(output_path, "cellphone_db_metadata.tsv"))

GetAssayData(tmp, "data") %>%
  DropletUtils::write10xCounts(
    path = here(output_path, "cellphone_db_data"),
    overwrite = TRUE
  )

file.rename(
  here(output_path, "cellphone_db_data/genes.tsv"),
  here(output_path, "cellphone_db_data/features.tsv")
)

rm(tmp)

dir.create(here(output_path, "cellphonedb_results"))

sys_cmd <- paste(
  here("lib/cellphonedb.sh"),
  here(output_path, "cellphonedb_results"),
  here(plots_path),
  argv$cores
)
debug(logger, paste0("The system command is: ", sys_cmd))
system(sys_cmd)

pdf(file = here(plots_path, "cellphone_plot.pdf"), height = 12, width = 18)
CellPhoneDotPlot(output_path = here(output_path, "cellphonedb_results"))
graphics.off()

cache_rds(
  expr = {
    cpdb <- ReadCellPhone(output_path = here(output_path, "cellphonedb_results")) %>%
      pivot_wider(names_from = clusters, values_from = mean) %>%
      select(-pvalue) %>%
      remove_missing(., vars = names(.)[2:length(.)]) %>%
      column_to_rownames(var = "pair") %>%
      CreateSeuratObject()

    cpdb <- ScaleData(cpdb)

    pc <- 100
    cpdb <- RunPCA(cpdb, features = rownames(cpdb), npcs = pc)
    cpdb <- JackStraw(cpdb, dims = pc)
    cpdb <- ScoreJackStraw(cpdb, dims = 1:pc)
    cpdb
  }, hash = file.info(here(output_path, "cellphonedb_results")),
  file = paste0(script_id, "cellphone_pca.rds")
)

JackStrawPlot(cpdb, dims = 1:(pc / 2)) + ElbowPlot(cpdb, ndims = (pc / 2)) & NoLegend()

good_pca_vars <- cpdb@reductions$pca@feature.loadings %>%
  as.data.frame() %>%
  abs() %>%
  slice_max(PC_1, n = 50) %>%
  rownames()
pca_vars <- cpdb@reductions$pca@feature.loadings %>%
  as.data.frame() %>%
  .[which(rownames(.) %in% good_pca_vars), ]

iplot <- plot_ly(
  data = pca_vars,
  x = ~PC_1,
  y = ~PC_2,
  text = rownames(pca_vars)
) %>%
  add_segments(x = 0, y = 0, xend = ~PC_1, yend = ~PC_2) %>%
  add_segments(x = 0, xend = 0, y = -0.1, yend = 0.1, color = I("grey"), line = list(dash = "dash")) %>%
  add_segments(y = 0, yend = 0, x = -0.1, xend = 0.1, color = I("grey"), line = list(dash = "dash")) %>%
  layout(
    shapes = list(
      type = "circle",
      xref = "x", x0 = -0.1, x1 = 0.1,
      yref = "y", y0 = -0.1, y1 = 0.1,
      color = I("grey"),
      line = list(dash = "dash")
    ),
    yaxis = list(scaleanchor = "x", zeroline = FALSE, showgrid = FALSE),
    xaxis = list(zeroline = FALSE, showgrid = FALSE),
    showlegend = FALSE,
    title = "CellPhoneDB Feature Loading Plot"
  ) %>%
  add_text(textposition = "top right")

htmlwidgets::saveWidget(as_widget(iplot), here(plots_path, "cellphone_feature_loading.html"))
if (rerun_cellphone) {
  seurat_tmp <- LoadH5Seurat(file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat"),
    assays = c("RNA")
  )
  DefaultAssay(seurat_tmp) <- "RNA"
  tmp <- subset(seurat_tmp, timepoint == "Diagnosis" & stemness == "Nonstem")
  tmp <- NormalizeData(tmp)

  Idents(diagnosis) %>%
    as_tibble(rownames = "Cell") %>%
    rename(cell_type = value) %>%
    write_tsv(here("data/cellphonedb_in/metadata.tsv"))

  GetAssayData(tmp, "data") %>%
    DropletUtils::write10xCounts(
      path = here("data/cellphonedb_in/data"),
      overwrite = TRUE
    )

  file.rename(
    here("data/cellphonedb_in/data/genes.tsv"),
    here("data/cellphonedb_in/data/features.tsv")
  )

  rm(tmp, seurat_tmp)

  system(here("lib/cellphonedb.sh"))
}
pdf(file = here(plots_path, "-cellphone_plot.pdf"), height = 12, width = 18)
CellPhoneDotPlot(output_path = here(output_path, "cellphonedb_results"))
graphics.off()

cache_rds(
  expr = {
    cpdb <- ReadCellPhone(output_path = here(output_path, "cellphonedb_results")) %>%
      pivot_wider(names_from = clusters, values_from = mean) %>%
      select(-pvalue) %>%
      remove_missing(., vars = names(.)[2:length(.)]) %>%
      column_to_rownames(var = "pair") %>%
      CreateSeuratObject()

    cpdb <- ScaleData(cpdb)

    pc <- 100
    cpdb <- RunPCA(cpdb, features = rownames(cpdb), npcs = pc)
    cpdb <- JackStraw(cpdb, dims = pc)
    cpdb <- ScoreJackStraw(cpdb, dims = 1:pc)
    cpdb
  }, hash = file.info(here(output_path, "cellphonedb_results")),
  file = paste0(script_id, "cellphone_pca.rds")
)

JackStrawPlot(cpdb, dims = 1:(pc / 2)) + ElbowPlot(cpdb, ndims = (pc / 2)) & NoLegend()

good_pca_vars <- cpdb@reductions$pca@feature.loadings %>%
  as.data.frame() %>%
  abs() %>%
  slice_max(PC_1, n = 50) %>%
  rownames()
pca_vars <- cpdb@reductions$pca@feature.loadings %>%
  as.data.frame() %>%
  .[which(rownames(.) %in% good_pca_vars), ]

iplot <- plot_ly(
  data = pca_vars,
  x = ~PC_1,
  y = ~PC_2,
  text = rownames(pca_vars)
) %>%
  add_segments(x = 0, y = 0, xend = ~PC_1, yend = ~PC_2) %>%
  add_segments(x = 0, xend = 0, y = -0.1, yend = 0.1, color = I("grey"), line = list(dash = "dash")) %>%
  add_segments(y = 0, yend = 0, x = -0.1, xend = 0.1, color = I("grey"), line = list(dash = "dash")) %>%
  layout(
    shapes = list(
      type = "circle",
      xref = "x", x0 = -0.1, x1 = 0.1,
      yref = "y", y0 = -0.1, y1 = 0.1,
      color = I("grey"),
      line = list(dash = "dash")
    ),
    yaxis = list(scaleanchor = "x", zeroline = FALSE, showgrid = FALSE),
    xaxis = list(zeroline = FALSE, showgrid = FALSE),
    showlegend = FALSE,
    title = "CellPhoneDB Feature Loading Plot"
  ) %>%
  add_text(textposition = "top right")

htmlwidgets::saveWidget(as_widget(iplot), here(plots_path, "cellphone_feature_loading.html"))
