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
parser$add_argument(
  "--batch-correct-method",
  "-M",
  "--bc-method",
  help = "which batch correction method to use",
  default = "CCA",
  choices = c("CCA", "RPCA")
)

argv <- parser$parse_args()

####### Define Functions #######################################################
DoDimensionReductions <- function(seurat,
                                  batch_vars = NULL,
                                  assay = NULL,
                                  clean = TRUE,
                                  scale = TRUE,
                                  features = NULL,
                                  resolution = 0.8) {
  assay <- assay %||% Seurat::DefaultAssay(seurat)
  if (clean) seurat <- Seurat::DietSeurat(seurat)

  verbosity <- function(x) {
    return(suppressMessages(x))
  }

  if (clean | scale) {
    info(logger, "Scaling data...")
    seurat <- verbosity(Seurat::FindVariableFeatures(seurat))
    features <- features %||% Seurat::VariableFeatures(seurat)

    if (Seurat:::IsSCT(seurat[[assay]])) {
      debug(logger, "SCTAssay object detected, getting residuals instead...")
      seurat <- verbosity(Seurat::GetResidual(seurat, features = features))
    } else {
      seurat <- verbosity(Seurat::ScaleData(seurat))
    }
  }
  info(logger, "Running PCA...")
  seurat <- Seurat::RunPCA(seurat, verbose = FALSE)
  reduction <- "pca"
  if (!is.null(batch_vars)) {
    info(logger, "Running Harmony...")
    seurat <- verbosity(harmony::RunHarmony(seurat, batch_vars, assay.use = assay))
    reduction <- "harmony"
  }
  info(logger, "Running UMAP...")
  seurat <- Seurat::RunUMAP(seurat, dims = 1:10, reduction = reduction, verbose = FALSE)

  info(logger, "Finding Clusters...")
  seurat <- verbosity(Seurat::FindNeighbors(seurat, dims = 1:19))
  seurat <- verbosity(Seurat::FindClusters(seurat, resolution = resolution))
  return(seurat)
}

FindClusterFreq <- function(object, metadata, cluster, sort_by = NULL) {
  # quote expressions
  md <- rlang::enquos(metadata)
  c <- rlang::enquo(cluster)
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

CalculatePValues.chi <- function(x) {
    if (is.matrix(x$obs)) {
        otmp <- apply(x$obs, 1, sum)
        etmp <- apply(x$exp, 1, sum)
    }
    else {
        otmp <- x$obs
        etmp <- x$exp
    }
    df <- (sum(1 * (etmp > 0))) - 1
    pval <- pchisq(x$chisq, df, lower.tail = FALSE)
    return(pval)
}

####### End Function Definitions ###############################################

if (argv$column_names) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratDisk)
    library(tidyverse)
  })
  seurat <- LoadH5Seurat(
    file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat"),
    assays = c("RNA")
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

# set plan
if (argv$cores == 0) argv$cores <- availableCores()[[1]]

if (argv$slurm) info(logger, "We're planning to use slurm for some futures")

SetPlan <- function(schedule = FALSE) {
  if (argv$slurm && schedule) {
    debug(logger, "The plan is set to batchtools using slurm")
    walltime_hours <- 24
    mem_gb <- 64
    resources <- list(memory = mem_gb, ncpus = 4, walltime = walltime_hours * 60^2)
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

furrr_options <- furrr_options(stdout = FALSE, seed = 1000000)

# so the plan here is to separate the computation into 4 steps 
# and use cache_rds to save the steps along the way
# Step 1: Subset object
# Step 2: If we want to integrate, Split data, and run SCTransform
# Step 3: If we want to integrate, integrate data using IntegrateData
# Step 4: Do dimension reductions

diagnosis <- cache_rds(
  expr = {
    assays_to_load <- c("RNA", "SCT")
    if (argv$batch_correct) assays_to_load <- c("RNA")
    assays_to_print <- glue_collapse(assays_to_load, sep = ", ", last = " and ")
    info(logger, glue("The assay(s) being loaded is/are {assays_to_print}"))
    seurat <- LoadH5Seurat(
      file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat"),
      assays = assays_to_load,
      verbose = FALSE
    )

    debug(logger, "Subsetting Seurat Object")
    MatchValuesForSubset <- function(val, col) {
      rownames(filter(seurat[[]], .data[[col]] == val))
    }
    debug(logger, glue("The class of cols is {class(argv$cols)}"))
    debug(logger, glue("The class of vals is {class(argv$vals)}"))
    res <- map(argv$cols, .f = function(col, vals) {
      unique(unlist(lapply(vals, MatchValuesForSubset, col)))
    }, argv$vals)

    cells <- reduce(res, .f = function(x, y) {
      x[match(x, y)]
    })
    cells <- cells[!is.na(cells)]

    diagnosis <- seurat[, cells]
    diagnosis
  },
  file = "seurat_subset.rds",
  dir = paste0(output_path, "/cache/"),
  rerun = argv$invalidate,
  hash = list(
    file.info(
      file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat")
    ),
    argv$vals,
    argv$cols
  )
)

if (argv$batch_correct) {
  object_list <- cache_rds(
    {
      debug(logger, "Splitting data by Patient ID")
      object_list <- SplitObject(diagnosis, split.by = "patient_id")
      info(logger, glue("Starting Seurat Integration on {length(object_list)} datasets"))
      debug(logger, "Normalizing using SCTransform")
      SetPlan(schedule = TRUE)
      object_list <- future_map(
        object_list,
        .f = SCTransform,
        .options = furrr_options,
        verbose = FALSE
      )
      object_list
    },
    file = "seurat_normalized.rds",
    dir = paste0(output_path, "/cache/"),
    rerun = argv$invalidate,
    hash = list(diagnosis)
  )

SetPlan()

bc_filename <- glue("seurat_integrated.rds")

info(logger, glue("Preparing to integrate using {argv$batch_correct_method"))

  diagnosis <- cache_rds(
    {
      debug(logger, "Selecting Integration Features")
      features <- SelectIntegrationFeatures(
        object.list = object_list,
        nfeatures = 2000,
        verbose = FALSE
      )
      debug(logger, "Preparing SCTransformed data for Integration")
      object_list <- PrepSCTIntegration(
        object.list = object_list,
        anchor.features = features,
        verbose = FALSE
      )
      if (argv$batch_correct_method == "RPCA") {
        debug(logger, "Running PCA on each dataset")
        object_list <- future_map(object_list, RunPCA, features = features)
      }
      bc_method <- if_else(argv$batch_correct_method == "CCA", "cca", "rpca")
      debug(logger, "Finding Anchors")
      anchors <- FindIntegrationAnchors(
        object.list = object_list,
        anchor.features = features,
        normalization.method = "SCT",
        reduction = bc_method,
        verbose = FALSE
      )
      info(logger, "Integrating Data")
      diagnosis <- IntegrateData(
        anchorset = anchors,
        normalization.method = "SCT",
        verbose = FALSE
      )
      diagnosis
    },
    file = bc_filename,
    dir = paste0(output_path, "/cache/"),
    rerun = argv$invalidate,
    hash = list(object_list)
  )
}

info(logger, "Doing Dimension Reductions")

diagnosis <- cache_rds(
  {
    batch_vars <- c("seq_batch", "sort_batch")
    if (argv$batch_correct) batch_vars <- NULL

    diagnosis <- DoDimensionReductions(
      diagnosis,
      batch_vars = batch_vars,
      resolution = argv$resolution,
      clean = !argv$batch_correct,
      scale = !argv$batch_correct
    )
    diagnosis[["clusters"]] <- Idents(diagnosis)
    diagnosis
  },
  file = glue("seurat_dimred.rds"),
  dir = paste0(output_path, "/cache/"),
  rerun = argv$invalidate,
  hash = list(diagnosis, argv$batch_correct)
)

info(logger, "Doing Frequency Analysis")
## Find DE Clusters ----
source(here("cli/lib/target_clinical.R"))

debug(logger, "Finding Cluster Frequency")

wilcox_clusters <- FindClusterFreq(
  diagnosis[[]],
  c("patient_id", "prognosis"),
  "clusters"
) %>%
  group_by(clusters) %>%
  summarise(
    pval = WilcoxTestSafe(freq ~ prognosis)$result$p.value, .groups = "keep"
  ) %>%
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

pdf(
  file = here(plots_path, "sc_cluster_survival_analysis.pdf"),
  width = 3,
  height = 5
)
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
