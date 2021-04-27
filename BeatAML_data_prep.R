beat_aml <- read_tsv(here("aml_ohsu_2018/data_RNA_Seq_expression_cpm.txt"))

rename(beat_aml, GeneSymbol = Hugo_Symbol) %>%
  mutate(Entrez_Gene_Id = NULL) %>%
  write_tsv(here("cibersort_in/beat_aml.txt"))
