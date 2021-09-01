#!/usr/bin/env Rscript
source("renv/activate.R")
library(argparse)

# parse args

parser <- ArgumentParser("Do Differential Clustering Analysis")
parser$add_argument(
  "--cols",
  "-c",
  help = "colums to subest",
  nargs = "+"
)
parser$add_argument(
  "--vals",
  "-a",
  help = "values columns should be. Can specify more than 1 value per column",
  nargs = "+"
)
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
  help = "ID to use for outputs",
  default = "time"
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
  "--column-names",
  "-C",
  help = "Just print column names of object",
  action = "store_true"
)
)

argv <- parse_args(parser)

if (argv$column_names) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratDisk)
    library(tidyverse)
  })
  seurat <- LoadH5Seurat(
    file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat"),
    assays = c("RNA", "SCT")
  )
  metadata <- seurat[[]] %>%
    select(where(is.character)) %>%
    lapply(unique) %>%
    lapply(paste, collapse = ", ")

  paste0(names(metadata), ": ", metadata) %>%
    paste(collapse = "\n") %>%
    cat()

  quit(save = "no")
}

if (!is.na(argv$cols) && !is.na(argv$subset)) {
  warning("Both columns and expression supplied, using expression only")
  argv$cols <- NA
  argv$vals <- NA
}

if (!is.na(argv$subset)) {
  stop("Do not use this argument: --subset")
}

if (!is.na(argv$cols) || !is.na(argv$vals)) {
  if (is.na(argv$cols) && is.na(argv$vals)) {
    stop("You must provide both the cols and vals argument")
  }
}
# load other packages

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(here)
  library(rsinglecell)
  library(readxl)
  library(biomaRt)
  library(survival)
  library(survminer)
  library(xfun)
  library(msigdbr)
  library(furrr)
  library(EnhancedVolcano)
  library(log4r)
  library(tidyverse)
})

source(here("lib/functions.R"))
logger <- logger(threshold = argv$verbose)

if (argv$id == "time") argv$id <- format(Sys.time(), format = "%m-%d-%Y[%H-%M]")

if (argv$dir == "") {
  output_path <- paste0("outs/", argv$id)
  plots_path <- paste0("plots/", argv$id)
} else {
  output_path <- paste0(parser$run_dir, "/outs/", argv$id)
  plots_path <- paste0(parser$run_dir, "/plots/", argv$id)
}

debug(logger, paste0("Writing Outputs to: ", here(output_path)))

dir.create(here(output_path), showWarnings = FALSE)
dir.create(here(plots_path), showWarnings = FALSE)

if (argv$invalidate) info("The cache will be invalidated")
diagnosis <- cache_rds(
  expr = {
    seurat <- LoadH5Seurat(
      file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat"),
      assays = c("RNA", "SCT")
    )

    DefaultAssay(seurat) <- "SCT"

    debug(logger, "Subsetting Seurat Object")
    if (is.na(argv$subset)) {
      MatchValuesForSubset <- function(val, col) {
        rownames(filter(seurat[[]], .data[[col]] == val))
      }

      res <- map(argv$cols, .f = function(col, vals) {
        unique(unlist(lapply(vals, MatchValuesForSubset, col)))
      }, argv$vals)

      cells <- reduce(res, .f = function(x, y) {
        x[match(x, y)]
      })
      cells <- cells[!is.na(cells)]

      diagnosis <- seurat[, cells]
    } else {
      diagnosis <- subset(seurat, parse(text = argv$subset))
    }

    info(logger, "Doing Dimension Reductions")
    diagnosis <- DoDimensionReductions(
      diagnosis,
      batch_vars = c("seq_batch", "sort_batch")
    )

    diagnosis[["clusters"]] <- Idents(diagnosis)
    diagnosis
  },
  file = "seurat.rds",
  dir = paste0(output_path, "/cache/"),
  rerun = argv$invalidate,
  hash = list(
    file.info(
      file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat")
    )
  )
)

info(logger, "Doing Frequency Analysis")
## Find DE Clusters ----
source(here("cli/lib/target_clinical.R"))

debug(logger, "Finding Cluster Frequency")

FindClusterFreq <- function(object, metadata, cluster, sort_by = NULL, .debug = FALSE) {
  # quote expressions
  md <- rlang::enquos(metadata)
  c <- rlang::enquo(cluster)
  if (.debug == TRUE) browser()
  # select meta.data columns from seurat object and calculate freq
  if (is.null(sort_by)) {
    freq <- object %>%
      dplyr::group_by(dplyr::across(.cols = all_of(c(metadata, cluster)))) %>%
      summarise(n = n(), .groups = "keep") %>%
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

WilcoxTestSafe <- purrr::safely(wilcox.test, otherwise = list(p.value = NA))

debug(logger, "Finding Cluster Frequency again")
wilcox_clusters <- FindClusterFreq(
  diagnosis[[]],
  c("patient_id", "prognosis"),
  "clusters"
) %>%
  group_by(clusters) %>%
  summarise(pval = WilcoxTestSafe(freq ~ prognosis)$result$p.value, .groups = "keep") %>%
  filter(pval <= argv$sig_level) %>%
  pull(clusters)

debug(logger, paste(wilcox_clusters, collapse = " "))

freq <- FindClusterFreq(
  diagnosis[[]],
  c("patient_id", "prognosis"),
  "clusters"
) %>%
  group_by(clusters) %>%
  select(patient_id, freq, clusters)

debug(logger, "Calculating Survival")
survival_models <- freq %>%
  right_join(clinical, by = c("patient_id" = "USI")) %>%
  mutate(
    status = if_else(`First Event` == "Relapse", 1, 0),
    cluster_risk = if_else(freq >= mean(freq), "High", "Low")
  ) %>%
  dplyr::select(
    `Event Free Survival Time in Days`,
    status,
    cluster_risk,
    patient_id,
    freq
  ) %>%
  #remove_missing() %>%
  filter(!near(freq, mean(freq))) %>%
  nest() %>%
  mutate(survival = map(
    data,
    ~ survdiff(
      Surv(`Event Free Survival Time in Days`, status) ~ cluster_risk,
      data = .x
    )
  )) %>%
  mutate(
    chi_sq = map_dbl(survival, ~ .x[["chisq"]]),
    p_val = map_dbl(survival, CalculatePValues.chi),
    log_p = -log(p_val),
    freq = map_dbl(data, ~ mean(.x[["freq"]])),
    wilcox = if_else(
      clusters %in% wilcox_clusters,
      "wilcox_sig",
      "not_wilcox_sig"
    ),
    chi_sig = if_else(p_val <= argv$sig_level, "chi_sig", "not_chi_sig"),
    sd = map_dbl(data, ~ sd(.x[["freq"]])),
    mean_diff = map_dbl(data, ReturnDifferences)
  ) %>%
  arrange(p_val)

saveRDS(survival_models, here(output_path, "cache/survival_models.rds"))

pdf(file = here(plots_path, "sc_cluster_survival_analysis.pdf"), width = 3, height = 5)
ggplot(data = survival_models, mapping = aes(x = mean_diff, y = log_p)) +
  geom_hline(yintercept = -log(argv$sig_level)) +
  geom_point() +
  theme_classic() +
  ggrepel::geom_label_repel(
    data = filter(survival_models, p_val <= argv$sig_level),
    mapping = aes(label = clusters)
  ) +
  xlab("Distance Between Poor (+) and Favorable (-)") +
  ylab("-log(P value)")
graphics.off()

source(here("lib/WriteInvocation.R"))
WriteInvocation(argv, output_path = paste0(output_path, "/invocation"))
