#!/usr/bin/env Rscript
source("renv/activate.R")
library(argparse)

# parse args

parser <- ArgumentParser("Do Differential Clustering Analysis")
parser$add_argument(
  "--sig-level",
  "-l",
  help = "p value cutoff",
  default = 0.05
)
parser$add_argument(
  "--dir",
  "-d",
  help = "path to run directory",
  default = ""
)
parser$add_argument(
  "--id",
  "-i",
  help = "ID to use for outputs"
)
parser$add_argument(
  "--invalidate",
  "-I",
  help = "should the xfun cache be invalidated?",
  action = "store_true"
)
parser$add_argument(
  "--verbose",
  "-v",
  help = "should messages be printed? One of: DEBUG, INFO, WARN, ERROR",
  default = "INFO"
)
parser$add_argument(
  "--resolution",
  "-r",
  help = "resolution to use",
  default = 0.8
)
parser$add_argument(
  "--batch-correct",
  "-X",
  help = "Should integration be done using Seurat",
  action = "store_true"
)
parser$add_argument(
  "--cores",
  "--threads",
  "-t",
  help = "number of processes to use. 0 (defualt) means all cores",
  default = 0,
  type = "integer"
)
parser$add_argument(
  "--slurm",
  "-S",
  help = "should batchtools be used with slurm",
  action = "store_true"
)

argv <- parser$parse_args()

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(here)
  suppressWarnings(library(rsinglecell))
  library(readxl)
  library(biomaRt)
  library(survival)
  library(survminer)
  library(xfun)
  library(msigdbr)
  library(future)
  library(furrr)
  library(future.batchtools)
  library(EnhancedVolcano)
  library(log4r)
  library(tidyverse)
  library(glue)
})

options(future.globals.maxSize = 32 * 1024^3)
source(here("lib/functions.R"))
logger <- logger(threshold = argv$verbose)

if (argv$dir == "") {
  output_path <- paste0("outs/", argv$id)
  plots_path <- paste0("plots/", argv$id)
} else {
  output_path <- paste0(parser$run_dir, "/outs/", argv$id)
  plots_path <- paste0(parser$run_dir, "/plots/", argv$id)
}

bc_ext <- ""
if (argv$batch_correct) bc_ext <- "bc_"

debug(logger, paste0("Writing Outputs to: ", here(output_path)))

if (argv$cores == 0) argv$cores <- availableCores()[[1]]

if (argv$slurm) info(logger, "We're planning to use slurm for some futures")

SetPlan <- function(schedule = FALSE, ncpu = 1, mem = 64) {
  if (argv$slurm && schedule) {
    info(logger, glue("The plan is set to batchtools using slurm with {ncpu} cpus."))
    walltime_hours <- 24
    mem_gb <- mem
    resources <- list(
      memory = mem_gb,
      ncpus = ncpu,
      walltime = walltime_hours * 60^2
    )
    plan(
      list(
        tweak("batchtools_slurm", resources = resources),
        "multicore"
      )
    )
  } else {
    debug(logger, "The plan is set to multisession")
    plan("multisession", workers = argv$cores)
    info(logger, glue("The number of workers is {nbrOfWorkers()[[1]]}"))
  }
}

if (argv$batch_correct) {
data_filename <- list.files(
  path = here(output_path, "cache/"),
  full.names = TRUE,
  pattern = "^seurat_integrated_bc_"
)
} else {
  data_filename <- list.files(
  path = here(output_path, "cache/"),
  full.names = TRUE,
  pattern = "^seurat_dimred_"
)
}

if (length(data_filename) > 1) {
  fatal(logger, "There is more than one cached object")
  quit(status = 1)
}

seurat <- readRDS(data_filename)
info(logger, "Splitting Object by prognosis groups")
object_list <- SplitObject(seurat, split.by = "prognosis")
object_list_names <- names(object_list)

SetPlan(schedule = TRUE)
furrr_options <- furrr_options(stdout = FALSE, seed = 1000000)

DoDimensionReductions <- function(object) {
  debug(logger, "Finding Variable Features")
  seurat <- FindVariableFeatures(object)
  debug(logger, "Running PCA")
  seurat <- RunPCA(seurat, npcs = 30, verbose = FALSE)
  debug(logger, "Running UMAP")
  seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:30, verbose = FALSE)
  debug(logger, "Finding Neighbors")
  seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:30)
  debug(logger, "Finding Clusters")
  seurat <- FindClusters(seurat, resolution = argv$resolution)
  return(seurat)
}

info(logger, "Starting Dimension Reduction and Clustering")
object_list <- cache_rds(
  future_map(object_list, DoDimensionReductions, .options = furrr_options),
  file = paste0(bc_ext, "dim_red_custom_de.rds"),
  rerun = argv$invalidate,
  dir = paste0(output_path, "/cache/")
)

info(logger, "Annotating Cluster Names")
cluster_ids <- map(object_list, Idents)
cluster_ids_ann <- map2(object_list_names, cluster_ids, ~ paste0(.x, "_", .y))

object_list <- map2(cluster_ids_ann, object_list, function(idents, object) {
  Idents(object) <- idents
  return(object)
})

seurat_merged <- reduce(object_list, merge)
VariableFeatures(seurat_merged) <- rownames(seurat_merged)

debug(logger, "Running PCA on Merged Data")
seurat_merged <- cache_rds(
  RunPCA(
    seurat_merged,
    npcs = 30,
    verbose = FALSE,
    reduction.name = "pca_merged"
  ),
  file = paste0(bc_ext, "pca_de_tests.rds"),
  rerun = argv$invalidate,
  dir = paste0(output_path, "/cache/")
)
debug(logger, "Running UMAP on Merged Data")
seurat_merged <- cache_rds(
  RunUMAP(
    seurat_merged,
    reduction = "pca_merged",
    dims = 1:30,
    reduction.name = "umap_merged",
    verbose = FALSE
  ),
  file = paste0(bc_ext, "umap_de_tests.rds"),
  rerun = argv$invalidate,
  dir = paste0(output_path, "/cache/")
)

pdf(here(plots_path, "UMAP_clusters_split.pdf"))
map(object_list, DimPlot, reduction = "umap")
DimPlot(seurat_merged, reduction = "umap_merged")
graphics.off()

info(logger, "Starting DE Testing")
SetPlan(schedule = TRUE, ncpu = 8)
debug(logger, "Doing DE Tests on seperate datasets")
markers_sep <- cache_rds(
  future_map(object_list, FindAllMarkers, method = "MAST"),
  file = paste0(bc_ext, "sep_de_contrasts.rds"),
  rerun = argv$invalidate,
  dir = paste0(output_path, "/cache/")
)
walk2(
  markers_sep,
  object_list_names,
  ~ write_tsv(.x, glue("{output_path}/{bc_ext}cluster_split_DE_{.y}.tsv"))
)

cluster_ids_ann_unique <- map(cluster_ids_ann, unique)
names(cluster_ids_ann_unique) <- object_list_names

results <- list()
for (group in object_list_names) {
  debug(logger, glue("Doing DE Tests on {group} comparing to ref"))
  group_sep <- Idents(seurat_merged) %>%
    str_subset(group) %>%
    unique()
  ref_cells <- unique(Idents(seurat_merged)[
    !(Idents(seurat_merged) %in% group_sep)
  ])
  debug(logger, "    Finished getting cells. Starting DE Testing")
  results[[group]] <- cache_rds(
    future_map(
      group_sep,
      ~ FindMarkers(
        seurat_merged,
        ident.1 = .x,
        ident.2 = ref_cells,
        method = "MAST",
        assay = "RNA"
      )
    ),
    file = glue("{bc_ext}ref_DE_tests_{group}.rds"),
    rerun = argv$invalidate,
    dir = paste0(output_path, "/cache/")
  )
  debug(logger, "    DE testing done")
  results[[group]] <- map2(
    results[[group]],
    cluster_ids_ann_unique[[group]],
    ~ mutate(.x, cluster = .y)
  ) %>%
    map(rownames_to_column, var = "gene")

  results[[group]] <- reduce(results[[group]], bind_rows)
}

info(logger, "Doing broader DE Tests")
broad_results <- cache_rds(
  FindMarkers(
    seurat_merged,
    ident.1 = "Poor",
    ident.2 = "Favorable",
    method = "MAST",
    assay = "RNA",
    group.by = "prognosis"
  ),
  file = glue("{bc_ext}ref_DE_tests_broad.rds"),
  rerun = argv$invalidate,
  dir = paste0(output_path, "/cache/")
)
broad_results <- mutate(broad_results, ref = "PvF", cluster = "PvF") %>% rownames_to_column(var = "gene")

debug(logger, "Cleaning up results.")
results <- map2(
  results,
  rev(names(cluster_ids_ann_unique)),
  ~ mutate(.x, ref = .y)
)
results[["PvF"]] <- broad_results
results <- reduce(results, bind_rows)

write_tsv(results, here(output_path, glue("{bc_ext}cluster_split_DE_master.tsv")))

source(here("lib/WriteInvocation.R"))
WriteInvocation(argv, output_path = here(output_path, "invocation"))
