library(Seurat)
library(tidyverse)

seurat <- readRDS("outs/flt3_cebpa_stem/cache/seurat_365ef0e0d8dbade1c5eb33690b90ab6a.rds")
source("cli/lib/prepare_clin_info.R")

markers <- FindAllMarkers(seurat, only.pos = TRUE, test.use = "LR")
colnames(markers)
markers %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice(1:25) %>%
  write_tsv("data/lsc_cluster_sig.tsv")