#!/usr/bin/env Rscript
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(glue)
library(progressr)

targeted_sequencing_dir <- "/fs/ess/PCCF0022/Datasets/AML_scRNAseq_dnw/targeted_sequencing/count_out"

libraries <- list.dirs(
    targeted_sequencing_dir,
    full.names = TRUE,
    recursive = FALSE
    )
libraries <- glue("{libraries}/outs/raw_feature_bc_matrix")

aml <- libraries[str_detect(libraries, "_aml")]
pan <- libraries[str_detect(libraries, "_pan_cancer")]
ids <- str_extract(pan, "[:alnum:]+_") %>% str_remove("_pan_cancer")

seurat <- LoadH5Seurat(
      "data/preprocessed.h5Seurat",
      verbose = FALSE
    )

file.copy("data/preprocessed.h5Seurat", "data/preprocessed.h5Seurat.bak")

split <- SplitObject(seurat, split.by = "library_id")
panels <- libraries

ids <- str_extract(panels, "/([:alnum:]|_)+_(aml|pan)") %>%
    str_remove("_(aml|pan)") %>%
    str_remove("/")

names <- str_extract(panels, "(aml|pan_cancer)")

params <- data.frame(
    panels = panels,
    ids = ids,
    names = names
)
panel_names <- unique(names)

a <- list()
handlers(global = TRUE)

aml_genes <- read_csv(
    "data/target_panel/AML_IDT.target_panel.csv",
    skip = 5,
    col_types = "__c"
    ) %>%
    separate(bait_id, into = c(NA, "Gene", NA)) %>%
    pull(Gene) %>%
    unique()

pan_cancer_genes <- read_csv(
    "data/target_panel/pan_cancer_v1.0_GRCh38-2020-A.target_panel.csv",
    skip = 5,
    col_types = "__c"
    ) %>%
    separate(bait_id, into = c(NA, "Gene", NA)) %>%
    pull(Gene) %>%
    unique()

PrepareTargetedCounts <- function(name, params, split) {
    iter_params <- params[params$names == name, ]
    panl <- iter_params$panels
    idx <- iter_params$ids
    p <- progressr::progressor(along = idx)
    f <- list()
    p(message = glue::glue("Preparing panel: {name}"), class = "sticky")
    for (i in seq_along(idx)) {
        s <- split[[idx[[i]]]]
        # just in case s does not exist
        if (!is.null(s)) {
            barcodes <- SeuratObject::Cells(s)
            striped_barcodes <- stringr::str_remove(barcodes, "-1_[:alnum:]+")
            suffix <- stringr::str_extract(barcodes, "-1_[:alnum:]+$")[[1]]
            p(amount = 0, message = glue::glue("Loading count matrix for {idx[[i]]}"))
            counts <- Seurat::Read10X(panl[[i]], strip.suffix = TRUE)
            filtered <- counts[ , striped_barcodes]
            colnames(filtered) <- glue::glue("{colnames(filtered)}{suffix}")
            f[[i]] <- filtered
        }
        p()
    }
    message(glue::glue("Preparing Assay Object"))
    merged_counts <- purrr::reduce(f, cbind)

    if (name == "aml") {
        merged_counts <- merged_counts[
            rownames(merged_counts) %in% aml_genes,
            ]
        }
    if (name == "pan_cancer") {
        merged_counts <- merged_counts[
            rownames(merged_counts) %in% pan_cancer_genes,
            ]
        }
    if (!(name %in% c("aml", "pan_cancer"))) {
        merged_counts <- merged_counts[rowSums(merged_counts) != 0, ]
        }
    assay <- SeuratObject::CreateAssayObject(counts = merged_counts)
    return(assay)
}

assays <- map(panel_names, PrepareTargetedCounts, params, split)
names(assays) <- stringr::str_to_upper(panel_names)

for (i in seq_along(assays)) {
    seurat[[names(assays)[[i]]]] <- assays[[i]]
}

seurat <- map(names(assays), ~ NormalizeData(seurat, assay = .x)))

seurat <- map(
    names(assays),
    ~ SCTransform(
        seurat,
        assay = .x,
        vars.to.regress = "percent.mt",
        new.assay.name = glue("{.x}_SCT")
        )
    )

SaveH5Seurat(
    seurat,
    "data/preprocessed.h5Seurat",
    overwrite = TRUE,
    verbose = FALSE
    )
