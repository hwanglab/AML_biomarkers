#! /bin/env Rscript
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

argv <- parser$parse_args()

# load more libraries
suppressPackageStartupMessages({
  library(here)
  library(log4r)
  library(tidyverse)
  library(glue)
})
logger <- logger(threshold = argv$verbose)
source(here("cli/lib/utils.R"))
output_path <- PrepareOutDir(argv)
StopIfOutputDirNotExist(output_path)

plots_path <- here(output_path, "plots")

decon <- readRDS(here(output_path, glue("cache/clinical_deconvoluted.rds")))


# decon$TCGA <- mutate(decon$TCGA, time = days_to_death)
decon$BeatAML <- mutate(decon$BeatAML, time = OS_DAYS)

target <- decon[c("TRAIN", "FLT3", "NEG", "CEBPA")]

TARGETPlot <- . %>%
  distinct(USI, .keep_all = TRUE) %>%
  pivot_longer(cols = starts_with("cluster"), names_to = "cluster", values_to = "frequency") %>%
  mutate(cluster = str_remove(cluster, "cluster") %>% as.numeric() %>% as.factor()) %>%
  ggplot(mapping = aes(x = USI, y = frequency, fill = cluster)) +
  geom_col() +
  theme_classic()

decon$BeatAML %>%
  pivot_longer(cols = starts_with("cluster"), names_to = "cluster", values_to = "frequency") %>%
  mutate(cluster = str_remove(cluster, "cluster") %>% as.numeric() %>% as.factor()) %>%
  ggplot(mapping = aes(x = Mixture, y = frequency, fill = cluster)) +
  geom_col() +
  theme_classic()
