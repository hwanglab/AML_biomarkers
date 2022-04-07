# function to turn a msigdb tibble into a list of gene sets with entrez gene ids
# this is the output from msigdb()

#' Turn Gene Set data.frame into a list
#'
#' @param df gene sets from `msigdb()`
#' @param gene_set should we retreive the Gene Set ID or the Gene Set Name?
#' @param gene_name should we use the Entrez IDs or Gene Symbols for gene names?
#' @export
#'
ListGeneSets <- function(df,
                         gene_set = c("id", "name"),
                         gene_name = c("entrez", "symbol")) {
  # set correct names
  if (gene_set == "id") {
    gene_set <- rlang::expr(gs_id)
  } else {
    gene_set <- rlang::expr(gs_name)
  }
  if (gene_name == "entrez") {
    gene_name <- rlang::expr(entrez_gene)
  } else {
    gene_name <- rlang::expr(gene_symbol)
  }
  # select the columns of interest
  gene_signatures <- dplyr::select(df, !!gene_set, !!gene_name)
  # nest tibble by signature
  list <- dplyr::nest_by(gene_signatures, !!gene_set)
  # set names (because that's what I did in my test)
  names <- list[[1]]
  names(list[[2]]) <- names
  # list the contents of the list
  data <- purrr::map(list[[2]], c, use.names = TRUE)
  # now unlist the list (will make it a single list, not a nested list)
  data <- unlist(data, recursive = FALSE)
  # set names
  names(data) <- names
  return(data)
}

RunGSVA <- function(seurat, gene_sets, assay = NULL, slot = "data", features = NULL,
                    average = TRUE, replicates = NULL) {
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
      seurat[["all_cells"]] <- "yes"
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

      all_cells <- AverageExpression(
        seurat,
        assay = assay,
        slot = slot,
        features = features,
        group.by = "all_cells"
      )[[1]]
      expr <- cbind(expr, all_cells)
      res <- GSVA::gsva(expr, gene_sets)
    } else {
      message("Averaging Expression in provided object...")
      expr <- Seurat::AverageExpression(seurat,
        assays = assay,
        slot = slot, features = features
      )[[1]]
      res <- GSVA::gsva(expr, gene_sets)
    }
  } else {
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
    res <- res[ , colnames(res) != "all"]
    all <- res[ , colnames(res) == "all"]
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
        top[[i]] <- limma::topTable(fit, n = Inf) %>%
          mutate(cluster = clusters[[i]])
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
        purrr::reduce(bind_rows) %>%
        mutate(AveExpr = 2^AveExpr)
    fake_top <- all %>%
      as_tibble() %>%
        set_names(c("AveExpr")) %>%
        rownames_to_column(var = "geneset") %>%
        mutate(
            logFC = 1,
            t = logFC, 
            P.Value = 0,
            adj.P.value = P.Value,
            B = AveExpr,
            cluster = "all"
        )
    top <- bind_rows(top, fake_top)
    plots <- patchwork::wrap_plots(plot)
    return(list(`Top Table` = top, plots = plot))
}