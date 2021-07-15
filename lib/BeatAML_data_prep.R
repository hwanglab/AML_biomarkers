library(here)
library(tidyverse)
library(rsinglecell)

beat_aml <- read_tsv(here("data/aml_ohsu_2018/data_RNA_Seq_expression_cpm.txt"))

a <- rename(beat_aml, GeneSymbol = Hugo_Symbol) %>%
  mutate(Entrez_Gene_Id = NULL) %>%
  transposeDF() %>%
  mutate(groups = ntile(x = row_number(), 5), .before = 1) %>%
  group_by(groups) %>%
  nest() %>%
  mutate(data = map(data, transposeDF, rowname_col = "GeneSymbol", colname_col = "samples"))


for (group in a$groups) {
  dat <- a$data[[group]]
  fname <- paste0("cibersort_in/beat_aml_", group, ".txt")
  write_tsv(dat, here(fname))
}
