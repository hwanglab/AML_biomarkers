library(SeuratDisk)
library(Seurat)

seurat <- LoadH5Seurat(
    "data/preprocessed.h5Seurat",
    assays = c("SCT", "RNA"),
    reductions = FALSE,
    graphs = FALSE,
    neighbors = FALSE,
    verbose = FALSE
)

genes <- c("CD117", "CD38", "CD34")
seurat <- GetResidual(seurat, genes)
VlnPlot(seurat, features = genes, group.by = "stemness", slot = "scale.data", pt.size = 0)

FindMarkers(seurat, features = c("CD38", "CD34"), group.by = "stemness", ident.1 = "Nonstem")

##      p_val avg_log2FC pct.1 pct.2 p_val_adj
## CD38     0  0.2508844 0.322 0.204         0
## CD34     0 -0.2541207 0.505 0.684         0
