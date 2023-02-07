library(Seurat)
library(here)
library(tidyverse)
library(readxl)
library(furrr)
library(patchwork)

seurat <- readRDS(here("outs/01_seurat_from_10x.rds"))

# Normalizing HTO
DefaultAssay(seurat) <- "HTO"
seurat <- NormalizeData(seurat, normalization.method = "CLR")

## Demux HTO ----
seurat <- HTODemux(seurat,
                   kfunc = "clara",
                   positive.quantile = 0.99,
                   nsamples = 500)

# remove negatives and doublets
table(seurat$HTO_classification.global, seurat$library_id) %>%
  t() %>%
  as.data.frame() %>%
  as_tibble() %>%
  rename(library_id = Var1, status = Var2) %>%
  pivot_wider(names_from = status, values_from = Freq) %>%
  rowwise() %>%
  mutate(
    total_counts = sum(Doublet, Negative, Singlet)
  ) %>%
  mutate(
    across(.cols = Doublet:Singlet, .fns = ~ .x / total_counts, .names = "{.col}_freq")
  ) %>%
  write_tsv("../biomarkers/one_off_scripts/per_library_stats.tsv")

