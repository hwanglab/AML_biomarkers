tcga_files <- list.files(here("tcga"), recursive = TRUE, pattern = ".gz$", full.names = TRUE)

tcga_data <- map(tcga_files, function(file_name) {
  folder <- here("tcga/")
  sam <- str_remove(file_name, folder)
  data <- read.table(file_name)
  colnames(data) <- c("ensembl_gene_id_version", sam)
  return(data)
})

tcga_data <- reduce(tcga_data, full_join)
tcga_clinical <- read_tsv(here("tcga/clinical.cart.2021-04-22/clinical.tsv"))
tcga_bridge <- read_tsv(here("tcga/gdc_sample_sheet.2021-04-22.tsv"))

tcga_ann <- inner_join(tcga_bridge, tcga_clinical, by = c("Case ID" = "case_submitter_id")) %>%
  unite(data_index, `File ID`, `File Name`, sep = "/")
tcga_ann_sam <- select(tcga_ann, data_index, `Case ID`)
tcga_sam <- colnames(tcga_data) %>%
  as.data.frame() %>%
  set_names("data_index") %>%
  left_join(tcga_ann_sam) %>%
  remove_missing() %>%
  distinct() %>%
  pull(`Case ID`) 
colnames(tcga_data) <- c("ensembl_gene_id_version", tcga_sam)

tcga_data <- EnsemblToHGNC(tcga_data, "ensembl_gene_id_version")
write_tsv(tcga_data, file = here("tcga_cibersort.tsv"))
