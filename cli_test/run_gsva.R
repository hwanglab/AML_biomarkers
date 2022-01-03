#!/usr/bin/env Rscript
source("renv/activate.R")
library(argparse)

# parse args
parser <- ArgumentParser("Do GSVA")
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
  "--replicates",
  "-n",
  help = "number of psudoreplicates to use during GSVA",
  default = 3
)
parser$add_argument(
  "--cores",
  "-c",
  help = "number of cores to use",
  default = 5
)
parser$add_argument(
  "--verbose",
  "-v",
  help = "should messages be printed? One of: DEBUG, INFO, WARN, ERROR",
  default = "INFO"
)

argv <- parser$parse_args()

# load more libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(here)
  suppressWarnings(
    library(rsinglecell)
  )
  library(msigdbr)
  library(furrr)
  library(EnhancedVolcano)
  library(tidyverse)
  library(log4r)
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
diagnosis <- readRDS(data_filename)

# set plan
if (argv$cores == 1) {
  plan("sequential")
} else {
  plan("multisession", workers = argv$cores)
}

gene_sets <- list(
  HALLMARK = msigdbr(species = "Homo sapiens", category = "H"),
  GO = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP"),
  REACTOME = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME"),
  KEGG = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG"),
  ONCO = msigdbr(species = "Homo sapiens", category = "C6")
) %>%
  map(~ ListGeneSets(.x, gene_set = "name", gene_name = "symbol"))

options("future.globals.maxSize" = 10 * 1024^3)

info(logger, "Running GSVA")

furr_options <- furrr_options(
  stdout = FALSE,
  seed = 1000000,
  conditions = structure("condition", exclude = c("message", "warning"))
)

debug(logger, paste0("The Default Assay is: ", DefaultAssay(diagnosis)))

RunGSVA <- function(seurat, gene_sets, assay = NULL, slot = "data", features = NULL,
                    average = TRUE, replicates = NULL) {
  PackageCheck <- rsinglecell:::PackageCheck
  if (!PackageCheck("GSVA", error = FALSE)) {
    stop("Please install GSVA.")
  }
  if (!is.list(gene_sets)) {
    stop(
      "Gene sets must be of type `list` not ",
      class(gene_sets)
    )
  }
  assay <- assay %||% Seurat::DefaultAssay(seurat)
  feautres <- features %||% unique(unlist(gene_sets))
  if (!is.null(replicates)) average <- TRUE
  if (is.null(replicates)) replicates <- 1
  if (average) {
    if (replicates > 1) {
      if (!PackageCheck("scorecard", error = FALSE)) {
        stop("Please install scorecard to do replicates")
      }
      DefaultAssay(seurat) <- assay

      # remove extra information from Seurat object to save space
      seurat <- Seurat::DietSeurat(seurat,
        counts = FALSE,
        assays = assay,
        dimreducs = FALSE,
        graphs = FALSE
      )
      rep_names <- paste0("rep", 1:replicates)
      idents <- as.data.frame(Idents(seurat))
      idents$idents <- rownames(idents)
      split_rats <- rep_len(1 / replicates, replicates)

      rlang::inform(paste0("Splitting Object to create ", replicates, " replicates."))
      meta_sub <- scorecard::split_df(idents, "Idents(seurat)",
        ratios = split_rats,
        name_dfs = rep_names
      )

      # Get Cell IDs for each object and subset the object
      cells <- purrr::map(meta_sub, ~ .x[["idents"]])
      objects <- purrr::map(cells, ~ subset(seurat, cells = .x))

      rlang::inform(paste0("Averaging expression for ", replicates, " replicates."))
      expr_list <- purrr::map(objects, ~ Seurat::AverageExpression(.x,
        assay = assay,
        slot = slot,
        features = features
      )[[1]])
      new_cols <- purrr::map2(expr_list, rep_names, ~ .x %>%
        as.data.frame() %>%
        colnames() %>%
        paste0("_", .y))
      expr_list <- purrr::map2(expr_list, new_cols, ~ `colnames<-`(.x, .y))

      rlang::inform(paste0("Preparing to run GSVA on ", replicates, " psudobulk replicates."))
      suppressMessages({
        expr <- expr_list %>%
          purrr::map(as.data.frame) %>%
          purrr::map(tibble::rownames_to_column) %>%
          purrr::reduce(dplyr::full_join) %>%
          column_to_rownames() %>%
          as.matrix()
      })

      res <- GSVA::gsva(expr, gene_sets)
    } else {
      message("Averaging Expression in provided object...")
      expr <- Seurat::AverageExpression(seurat,
        assays = assay,
        slot = slot, features = features
      )[[1]]
      res <- GSVA::gsva(expr, gene_sets)
    }
  }
  else {
    stop("GSVA on single-cell expression has not been implemented well. Please run with average = TRUE only and make use of replicates")
    message("Getting Assay data for provided features...")
    dat <- Seurat::GetAssayData(seurat, assay = assay, slot = slot)
    res <- list(0)
    pb <- progress::progress_bar$new(length(gene_sets),
      format = "Running GSVA [:bar] :percent (:elapsed)"
    )
    for (i in seq_along(gene_sets)) {
      expr <- dat[which(dat@Dimnames[[1]] %in% unique(unlist(gene_sets[[i]]))), ]
      res[[i]] <- GSVA::gsva(as.matrix(expr), gene_sets[i])
      pb$tick()
    }
  }
  return(res)
}

RunStats <- function(res, pathway, p = 0.05, fc = 0.25) {
  clusters <- colnames(res) %>%
    tibble::as_tibble_col(column_name = "names") %>%
    tidyr::separate(names, into = c("cluster", "rep"), sep = "_") %>%
    dplyr::pull(cluster) %>%
    unique()
  top <- list(0)
  plot <- list(0)
  for (i in seq_along(clusters)) {
    new_names <- colnames(res) %>%
      tibble::as_tibble_col(column_name = "names") %>%
      tidyr::separate(names, into = c("cluster", "rep"), sep = "_") %>%
      dplyr::group_by(
        cluster2 = ifelse(
          cluster == clusters[[i]],
          yes = paste0("cluster", clusters[[i]]),
          no = paste0("not", clusters[[i]])
        )
      ) %>%
      dplyr::mutate(rep2 = dplyr::row_number())
    design <- model.matrix(~ 0 + new_names$cluster2)
    colnames(design) <- colnames(design) %>% str_remove("new_names\\$cluster2")
    con_str <- c(paste0(
      "cluster", clusters[[i]], " - not",
      clusters[[i]]
    ))
    contrasts <- limma::makeContrasts(
      contrasts = con_str,
      levels = design
    )
    fit <- limma::lmFit(res, design)
    fit <- limma::contrasts.fit(fit, contrasts = contrasts)
    fit <- limma::eBayes(fit)
    top[[i]] <- limma::topTable(fit, n = Inf) %>% mutate(cluster = clusters[[i]])
    labs <- stringr::str_remove(rownames(top[[i]]), "^(GOBP|HALLMARK|KEGG|REACTOME)_")
    plot[[i]] <- EnhancedVolcano::EnhancedVolcano(top[[i]],
      lab = labs, x = "logFC", y = "adj.P.Val", subtitle = paste0(
        "Cluster ",
        clusters[[i]]
      ), title = paste0(
        "GSVA results for ",
        pathway
      ), pCutoff = p, FCcutoff = fc, xlim = c(
        -1,
        1
      ), labSize = 1.75, drawConnectors = TRUE
    )
  }
  top <- top %>%
    purrr::map(rownames_to_column, var = "geneset") %>%
    purrr::reduce(bind_rows)
  plots <- patchwork::wrap_plots(plot)
  return(list(`Top Table` = top, plots = plot))
}
RunGSVAQuietly <- function(...) {
  f <- purrr::quietly(RunGSVA)
  res <- f(...)$result
  return(res)
}

gsva_res <- future_map(
  gene_sets,
  ~ RunGSVAQuietly(
    diagnosis,
    gene_sets = .x,
    replicates = argv$replicates
  ),
  .options = furr_options
)

info(logger, "GSVA Done")

RunStatsQuietly <- function(...) {
  f <- purrr::quietly(RunStats)
  res <- f(...)$result
  return(res)
}

stat_res <- future_map2(
  gsva_res,
  names(gsva_res),
  RunStatsQuietly,
  p = 0.05,
  .options = furr_options
)

file_names <- paste0(names(gene_sets), "_GSVA_volcano_plots.pdf")

top_tables <- map(stat_res, ~ .x[["Top Table"]])

top_tables %>%
  reduce(bind_rows) %>%
  list("All Gene Sets" = .) %>%
  append(top_tables) %>%
  openxlsx::write.xlsx(file = here(output_path, "GSVA_DE_results.xlsx"))

plots <- map(stat_res, ~ .x[["plots"]])

walk2(file_names, plots, .f = function(file, plot) {
  pdf(file = here(plots_path, file), onefile = TRUE)
  print(plot)
  graphics.off()
})

source(here("lib/WriteInvocation.R"))
WriteInvocation(argv, output_path = here(output_path, "invocation"))
