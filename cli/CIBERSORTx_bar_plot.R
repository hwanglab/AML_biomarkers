#! /bin/env Rscript

source("renv/activate.R")
library(argparse)

# parse args
parser <- ArgumentParser("Prepare Clinical Data")
parser$add_argument(
  "--dir",
  "-d",
  help = "path to run directory",
  default = ""
)
parser$add_argument(
  "--id",
  "-i",
  help = "ID to use for outputs"
)
parser$add_argument(
  "--verbose",
  "-v",
  help = "should messages be printed? One of: DEBUG, INFO, WARN, ERROR",
  default = "INFO"
)
parser$add_argument(
  "--test-results",
  "-t",
  help = "name of test for lasso model",
  required = TRUE
)

argv <- parser$parse_args()

# load more libraries
suppressPackageStartupMessages({
  library(here)
  library(log4r)
  library(tidyverse)
  library(glue)
  library(viridis)
  library(gghighlight)
})

if (argv$dir == "") {
  output_path <- paste0("outs/", argv$id)
  plots_path <- paste0("plots/", argv$id)
} else {
  output_path <- paste0(parser$run_dir, "/outs/", argv$id)
  plots_path <- paste0(parser$run_dir, "/plots/", argv$id)
}
if (!dir.exists(here(output_path))) {
  fatal(logger, "Output directory does not exist")
}

logger <- logger(threshold = argv$verbose)

bc_ext <- "no_bc"
if (argv$batch_correct) bc_ext <- argv$batch_correct_method

decon <- readRDS(here(output_path, glue("cache/clinical_deconvoluted.rds")))
target <- decon[c("FLT3", "NEG", "CEBPA")]

PrepareTARGETData <- . %>%
  distinct(USI, .keep_all = TRUE) %>%
  pivot_longer(cols = starts_with("cluster"), names_to = "cluster", values_to = "frequency") %>%
  mutate(cluster = str_remove(cluster, "cluster") %>% as.numeric() %>% as.factor())


FinishGGPlot <- function(highlight = NULL) {
  list(
    geom_col(),
    theme_classic(),
    scale_fill_viridis(discrete = TRUE),
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ),
    if (!is.null(highlight)) {
      gghighlight(cluster %in% highlight, calculate_per_facet = TRUE)
    } else {
      NULL
    }
  )
}

debug(logger, "Preparing Data")

target_data <- map(target, PrepareTARGETData)
beat_aml_data <- decon$BeatAML %>%
  pivot_longer(cols = starts_with("cluster"), names_to = "cluster", values_to = "frequency") %>%
  mutate(cluster = str_remove(cluster, "cluster") %>% as.numeric() %>% as.factor())

tcga_data <- decon$TCGA %>%
  pivot_longer(cols = starts_with("cluster"), names_to = "cluster", values_to = "frequency") %>%
  mutate(cluster = str_remove(cluster, "cluster") %>% as.numeric() %>% as.factor())


debug(logger, "Getting LASSO results")
lasso <- read_tsv(
  here(glue("{output_path}/{argv$test}/lasso_model_coefs.tsv")), 
  col_types = cols()
  )

lasso <- lasso[lasso$s0 != 0, ]
lasso <- lasso[2:nrow(lasso), ]

clusters <- lasso$rownames %>% str_remove("cluster") %>% as.numeric()
clusters <- clusters[!is.na(clusters)]

debug(logger, "Plotting Data")
plot <- map(target_data, ~ .x %>%
  ggplot(mapping = aes(x = USI, y = frequency, fill = cluster)) +
  facet_grid(~`First Event`, scales = "free_x", space = "free_x") +
  FinishGGPlot())

plot2 <- list(
  TCGA = ggplot(data = plot_tcga_data, mapping = aes(x = Mixture, y = frequency, fill = cluster)) +
    facet_grid(~vital_status, scales = "free_x", space = "free_x") +
    FinishGGPlot(),
  BeatAML = ggplot(data = plot_tcga_data, mapping = aes(x = Mixture, y = frequency, fill = cluster)) +
    facet_grid(~vital_status, scales = "free_x", space = "free_x") +
    FinishGGPlot()
)
plot <- splice(plot, plot2)
plot <- map2(plot, names(plot), ~ .x + labs(title = .y))

plot2 <- map(plot, ~ .x + FinishGGPlot(highlight = clusters))

pdf(glue("{plots_path}/CIBERSORTx_bar_plot.pdf"))
plot
graphics.off()

pdf(glue("{plots_path}/CIBERSORTx_bar_plot_highlight.pdf"))
plot2
graphics.off()
