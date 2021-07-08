
ReturnDifferences <- function(x) {
  x %>%
  #mutate(freq = scale(freq, center = FALSE)) %>%
  dplyr::group_by(status) %>%
  dplyr::summarise(freq_mean = mean(freq)) %>%
  dplyr::ungroup() %>%
  dplyr::summarise(value = max(freq_mean) - min(freq_mean),
            value2 = status[which.max(freq_mean)]) %>%
  dplyr::mutate(value3 = if_else(value2 == 1, value, -1 * value)) %>%
  dplyr::pull(value3)
}





CellPhoneDotPlot <- function(selected_rows = NULL,
                             selected_columns = NULL,
                             output_path = NULL,
                             means_separator = '\t',
                             pvalues_separator = '\t'
){
  if (is.null(output_path)) stop("Must provide output path")
  pvalues_path <- paste0(output_path, "/pvalues.txt")
  means_path <- paste0(output_path, "/means.txt")
  
  message("Reading Outputs...")
  all_pval <- read.table(pvalues_path, header = T, stringsAsFactors = F, sep = means_separator, comment.char = '', check.names = F)
  all_means <- read.table(means_path, header = T, stringsAsFactors = F, sep = pvalues_separator, comment.char = '', check.names = F)
  
  intr_pairs <- all_pval$interacting_pair
  all_pval <- all_pval[,-c(1:11)]
  all_means <- all_means[,-c(1:11)]
  
  if(is.null(selected_rows)){
    selected_rows <- intr_pairs
  }
  
  if(is.null(selected_columns)){
    selected_columns <- colnames(all_pval)
  }
  
  sel_pval <- all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means <- all_means[match(selected_rows, intr_pairs), selected_columns]
  
  message("Preparing Plot Data...")
  
  df_names <- expand.grid(selected_rows, selected_columns)
  pval <- unlist(sel_pval)
  pval[pval == 0] <- 0.0009
  plot.data <- cbind(df_names,pval)
  pr <- unlist(as.data.frame(sel_means))
  pr[pr == 0] <- 1
  plot.data <- cbind(plot.data,log2(pr))
  colnames(plot.data) <- c('pair', 'clusters', 'pvalue', 'mean')
  
  good_features <- plot.data %>%
    tidyr::separate(clusters, into = c("c1", "c2")) %>%
    dplyr::filter(c1 != c2 & pvalue <= 0.05 & abs(mean) > 3.5) %>%
    tidyr::unite(clusters, c1, c2, sep = "|") %>%
    dplyr::select(pair, clusters)
  
  plot.data <- plot.data[plot.data$pair %in% good_features$pair, ]
  plot.data <- plot.data[plot.data$clusters %in% good_features$clusters, ]
  
  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha = TRUE)(n = 399)
  
  message("Building Plot...")
  
  plot <- ggplot2::ggplot(plot.data, ggplot2::aes(x = clusters, y = pair)) +
    ggplot2::geom_point(ggplot2::aes(size = -log10(pvalue), color = mean)) +
    ggplot2::scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors = my_palette) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          axis.text = ggplot2::element_text(size = 14, colour = "black"),
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
          axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
          axis.title = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(size = 0.7, linetype = "solid", colour = "black"))
  
  return(plot)
}


ReadCellPhone <- function(selected_rows = NULL,
                          selected_columns = NULL,
                          output_path = NULL,
                          means_separator = '\t',
                          pvalues_separator = '\t'
){
  if (is.null(output_path)) stop("Must provide output path")
  pvalues_path <- paste0(output_path, "/pvalues.txt")
  means_path <- paste0(output_path, "/means.txt")
  
  all_pval <- read.table(pvalues_path, header = T, stringsAsFactors = F, sep = means_separator, comment.char = '', check.names = F)
  all_means <- read.table(means_path, header = T, stringsAsFactors = F, sep = pvalues_separator, comment.char = '', check.names = F)
  
  intr_pairs <- all_pval$interacting_pair
  all_pval <- all_pval[,-c(1:11)]
  all_means <- all_means[,-c(1:11)]
  
  if (is.null(selected_rows)){
    selected_rows <- intr_pairs
  }
  
  if (is.null(selected_columns)){
    selected_columns <- colnames(all_pval)
  }
  
  sel_pval <- all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means <- all_means[match(selected_rows, intr_pairs), selected_columns]
  
  df_names <- expand.grid(selected_rows, selected_columns)
  pval <- unlist(sel_pval)
  pval[pval == 0] <- 0.0009
  plot.data <- cbind(df_names,pval)
  pr <- unlist(as.data.frame(sel_means))
  pr[pr == 0] <- 1
  plot.data <- cbind(plot.data,log2(pr))
  colnames(plot.data) <- c('pair', 'clusters', 'pvalue', 'mean')
  
  return(plot.data)
}

PCAVariance <- function(seurat, assay = "RNA") {
  mat <- Seurat::GetAssayData(seurat, assay = assay, slot = "scale.data")
  pca <- seurat[["pca"]]
  
  # Get the total variance:
  total_variance <- sum(matrixStats::rowVars(mat))
  
  eigValues <- (pca@stdev)^2  ## EigenValues
  varExplained <- eigValues / total_variance
  return(varExplained)
}

#' Do a batch survival analysis
#' 
#' @param data named list of datasets to test on
#' @param time columns where survival times are stored
#' @param lasso_model a lasso model
#' @param event_col a named vector where the names are columns and the values are strings of desired outcomes
#' 
#' @return a list with plots and data.frame of statitics
#' 
#' @examples 
#' /dontrun{
#' # setup times in the data
#' times <- c("Event Free Survival Time in Days", "Overall Survival Time in Days", "days_to_death", "OS_DAYS")
#' 
#' # set up what the event colums are and what the events are labelled as
#' event_col <- c("First Event" = "relapse", "vital_status" = "Dead", "status" = 1)
#' 
#' # run the survival analysis
#' results <- BatchSurvival(data, times, lasso_model, event_col)
#' }
BatchSurvival <- function(data, time, lasso_model, event_col) {

    events_ <- rlang::syms(names(event_col))
    times <- c("Event Free Survival Time in Days", "Overall Survival Time in Days", "days_to_death", "OS_DAYS")

    stats_final <- list(0)
    plots_final <- list(0)
    for (dataset in seq_along(data)) {
        x <- as_tibble(data[[dataset]])
        avail_cols <- events_[names(event_col) %in% colnames(data[[dataset]])]
        avail_events <- event_col[names(event_col) %in% colnames(data[[dataset]])]

        stats <- list(0)
        plots <- list(0)
        
        for (st in seq_along(avail_cols)) {
            message("Using ", avail_cols[[st]], " on ", names(data)[[dataset]])
            x[["status"]] <- if_else(x[[avail_cols[[st]]]] == avail_events[[st]], 1, 0)
            x[["score"]] <- UseLASSOModelCoefs(x, coef(lasso_model))
            x[["score_bin"]] <- if_else(x$score >= median(x$score, na.rm = TRUE), "High", "Low")
            res <- list(0)

            for (time_axis in seq_along(times)) {
                pSuv <- purrr::possibly(DoSurvialAnalysis, otherwise = NA)
                desc <- paste0(names(data)[[dataset]], ": ", avail_cols[[st]])
                res[[time_axis]] <- pSuv(x,
                    !!rlang::sym(times[[time_axis]]),
                    status,
                    score,
                    group_by = score_bin,
                    description = desc,
                    lasso = lasso_model)
            }
            names(res) <- times
            res <- purrr::discard(res, ~ all(is.na(.x)))

            if (length(res) > 0) {
                stats_unclean <- map(res, ~ .x[["stat"]])
                plots[[st]] <- map(res, ~ .x[["plot"]])
                stats[[st]] <- stats_unclean %>%
                    map(broom::tidy, conf.int = TRUE) %>%
                    reduce(bind_rows) %>%
                    mutate(meta = names(res),
                        data = names(data)[[dataset]],
                        col = rep_len(rlang::as_name(avail_cols[[st]]), length(stats_unclean)))
            }
        }
        plots_final[[dataset]] <- plots
        stats_final[[dataset]] <- reduce(stats, bind_rows)
    }
    res1 <- stats_final %>% keep(is.data.frame) %>% reduce(bind_rows)    
    res2 <- plots_final %>% unlist(recursive = FALSE) %>% unlist(recursive = FALSE)
    return(list(stats = res1, plots = res2))
}
