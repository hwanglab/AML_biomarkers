
ReturnDifferences <- . %>%
  #mutate(freq = scale(freq, center = FALSE)) %>%
  group_by(status) %>%
  summarise(freq_mean = mean(freq)) %>%
  ungroup() %>%
  summarise(value = max(freq_mean) - min(freq_mean),
            value2 = status[which.max(freq_mean)]) %>%
  mutate(value3 = if_else(value2 == 1, value, -1 * value)) %>%
  pull(value3)





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
    separate(clusters, into = c("c1", "c2")) %>%
    filter(c1 != c2 & pvalue <= 0.05 & abs(mean) > 3.5) %>%
    unite(clusters, c1, c2, sep = "|") %>%
    select(pair, clusters)
  
  plot.data <- plot.data[plot.data$pair %in% good_features$pair, ]
  plot.data <- plot.data[plot.data$clusters %in% good_features$clusters, ]
  
  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha = TRUE)(n = 399)
  
  message("Building Plot...")
  
  plot <- ggplot(plot.data, aes(x = clusters, y = pair)) +
    geom_point(aes(size = -log10(pvalue), color = mean)) +
    scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors = my_palette) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_text(size = 14, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size = 12, colour = "black"),
          axis.title = element_blank(),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
  
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
