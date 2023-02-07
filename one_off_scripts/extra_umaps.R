
library(Seurat)
library(tidyverse)

# Load the data
seurat <- readRDS("outs/flt3_cebpa_stem/cache/seurat_365ef0e0d8dbade1c5eb33690b90ab6a.rds") 
all_cells <- readRDS("outs/all_cells/cache/seurat_dimred_ef84df7a0aceafed6e8d178dd89b190e.rds")
nonstem <- readRDS("outs/flt3_cebpa/cache/seurat_d307f778811553cd7be8a49668d726e4.rds")

Idents(seurat) <- "patient_id"
low_med_read_pats <- c("PAUJMC", "PASXVC")
names(low_med_read_pats) <- low_med_read_pats
low_med_read_cells <- map(
    low_med_read_pats, 
    ~ possibly(WhichCells, otherwise = NULL)(seurat, ident = .x)
)

high_doublet_libs <- c("PASDET", "PARXEC", "PARVEX", "PAUJMC", "PASXVC", "PAUJMC")
names(high_doublet_libs) <- high_doublet_libs
high_doublet_cells <- map(
    high_doublet_libs, 
    ~ possibly(WhichCells, otherwise = NULL)(seurat, ident = .x)
)

Idents(seurat) <- "patient_id"
high_doublet_libs <- c(
    "PAUSFM", "PAUXYG", "PATHZC", "PATIAB", "PAWTWW", "PAUNAZ", "PAXMKT", 
    "PAWAMB", "PARDZW", "PARSGS", "PARYCG", "PARZVA",  "PAUJMC", "PASXVC", 
    "PAUJMC"
)
names(high_doublet_libs) <- high_doublet_libs
high_doublet_cells_nonstem <- map(
    high_doublet_libs, 
    ~ possibly(WhichCells, otherwise = NULL)(nonstem, ident = .x)
)

pdf("one_off_scripts/extra_umaps.pdf", width = 10, height = 10)
print(UMAPPlot(seurat, cells.highlight = low_med_read_cells))
UMAPPlot(
    seurat, 
    cells.highlight = high_doublet_cells, 
    cols.highlight = c("#D30C7B", "#3A2D32", "#A2AD91"), 
    pt.size = 0.5
)
UMAPPlot(nonstem, cells.highlight = high_doublet_cells_nonstem)
UMAPPlot(all_cells, group.by = "stemness")
graphics.off()
