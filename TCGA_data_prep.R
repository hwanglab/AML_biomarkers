tcga_files <- list.files(here("data/tcga"), recursive = TRUE, pattern = ".gz$", full.names = TRUE)

tcga_data <- map(tcga_files, function(file_name) {
  folder <- here("data/tcga/")
  sam <- str_remove(file_name, folder)
  data <- read.table(file_name)
  colnames(data) <- c("ensembl_gene_id_version", sam)
  return(data)
})

tcga_data <- reduce(tcga_data, full_join)
tcga_clinical <- read_tsv(here("data/tcga/clinical.cart.2021-04-22/clinical.tsv"), na = "'--")
tcga_bridge <- read_tsv(here("data/tcga/gdc_sample_sheet.2021-04-22.tsv"))

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

tcga_ann2 <- read_tsv(here("clinical_info/nationwidechildrens.org_clinical_patient_laml.txt"), skip = 1, na = c("[Not Available]", "[Not Applicable]")) %>%
  filter(bcr_patient_uuid != "CDE_ID:") %>%
  map_df(parse_guess) %>%
  mutate(flt3_status = str_extract(molecular_analysis_abnormality_testing_result, "FLT3 Mutation [:alpha:]+")) %>%
  separate(flt3_status, into = c(NA, NA, "flt3_status"), sep = " ") %>%
  rename(`Case ID` = bcr_patient_barcode) %>%
  left_join(tcga_ann) %>%
  rename(case_submitter_id = `Case ID`)
