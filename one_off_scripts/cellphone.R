RunCellPhoneDBAnalysis <- function(output_path, script_id, idents, plots_path, cells = NULL) {
  
  `%>%` <- magrittr::`%>%`
  `%||%` <- SeuratObject::`%||%`
  source(here::here("lib/functions.R"))
  seurat_tmp <- SeuratDisk::LoadH5Seurat(file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat"),
    assays = c("RNA")
  )
  Seurat::DefaultAssay(seurat_tmp) <- "RNA"
  cells <- cells %||% colnames(seurat_tmp)
  tmp <- seurat_tmp[, cells]
  tmp <- Seurat::NormalizeData(tmp)

  idents %>%
    tibble::as_tibble(rownames = "Cell") %>%
    dplyr::rename(cell_type = value) %>%
    readr::write_tsv(here::here("data/cellphonedb_in/metadata.tsv"))

  Seurat::GetAssayData(tmp, "data") %>%
    DropletUtils::write10xCounts(
      path = here::here("data/cellphonedb_in/data", script_id),
      overwrite = TRUE
    )

  file.rename(
    here::here("data/cellphonedb_in/data", script_id, "genes.tsv"),
    here::here("data/cellphonedb_in/data", script_id, "features.tsv")
  )

  rm(tmp, seurat_tmp)

  system(here::here("lib/cellphonedb.sh"))

  pdf(file = here::here(plots_path, "-cellphone_plot.pdf"), height = 12, width = 18)
  CellPhoneDotPlot(output_path = here::here(output_path, "cellphonedb_results"))
  graphics.off()

  xfun::cache_rds(
    expr = {
      cpdb <- ReadCellPhone(output_path = here::here(output_path, "cellphonedb_results")) %>%
        tidyr::pivot_wider(names_from = clusters, values_from = mean) %>%
        dplyr::select(-pvalue) %>%
        ggplot2::remove_missing(., vars = names(.)[2:length(.)]) %>%
        tibble::column_to_rownames(var = "pair") %>%
        Seurat::CreateSeuratObject()

      cpdb <- Seurat::ScaleData(cpdb)

      pc <- 100
      cpdb <- Seurat::RunPCA(cpdb, features = rownames(cpdb), npcs = pc)
      cpdb <- Seurat::JackStraw(cpdb, dims = pc)
      cpdb <- Seurat::ScoreJackStraw(cpdb, dims = 1:pc)
      cpdb
    }, hash = file.info(here::here(output_path, "cellphonedb_results")),
    file = paste0(script_id, "cellphone_pca.rds")
  )
  # JackStrawPlot(cpdb, dims = 1:(pc / 2)) + ElbowPlot(cpdb, ndims = (pc / 2)) & NoLegend()

  good_pca_vars <- cpdb@reductions$pca@feature.loadings %>%
    as.data.frame() %>%
    abs() %>%
    dplyr::slice_max(PC_1, n = 50) %>%
    rownames()
  pca_vars <- cpdb@reductions$pca@feature.loadings %>%
    as.data.frame() %>%
    .[which(rownames(.) %in% good_pca_vars), ]

  iplot <- plotly::plot_ly(
    data = pca_vars,
    x = ~PC_1,
    y = ~PC_2,
    text = rownames(pca_vars)
  ) %>%
    plotly::add_segments(x = 0, y = 0, xend = ~PC_1, yend = ~PC_2) %>%
    plotly::add_segments(x = 0, xend = 0, y = -0.1, yend = 0.1, color = I("grey"), line = list(dash = "dash")) %>%
    plotly::add_segments(y = 0, yend = 0, x = -0.1, xend = 0.1, color = I("grey"), line = list(dash = "dash")) %>%
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
    plotly::add_text(textposition = "top right")

  htmlwidgets::saveWidget(as_widget(iplot), here(plots_path, "cellphone_feature_loading.html"))

  return(invisible())
}
