library(Seurat)
library(SeuratDisk)
library(here)
library(rsinglecell)
library(readxl)
library(biomaRt)
library(survival)
library(survminer)
library(furrr)
library(tidyverse)

seurat <- LoadH5Seurat(file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat"),
  assays = c("ADT", "SCT")
)

stem <- readRDS(file.path(Sys.getenv("AML_DATA"), "11-stem.rds"))

stem_paired <- stem[[]] %>% filter(paired == "T", patient_id != "PAUVEF")
FindClusterFreqDF(stem_paired, c("patient_id", "timepoint"), "clusters") %>%
  group_by(cluster, .drop = FALSE) %>%
  summarise(
    pval = wilcox.test(freq ~ timepoint)$p.value,
    .groups = "keep"
  ) %>%
  arrange(pval)

stem_diagnosis <- subset(stem, timepoint == "Diagnosis")
FindClusterFreq(stem_diagnosis, c("patient_id", "prognosis"), "clusters") %>%
  group_by(cluster) %>%
  summarise(pval = wilcox.test(freq ~ prognosis)$p.value, .groups = "keep") %>%
  arrange(pval)

test <- FindClusterFreq(stem, c("patient_id", "timepoint"), "clusters") %>%
  pivot_wider(names_from = cluster, values_from = freq, names_prefix = "cluster")

sapply(test[3:23], sd)

### New Biomarkers

diagnosis <- subset(seurat, timepoint == "Diagnosis")
diagnosis <- DoDimensionReductions(diagnosis, batch_vars = c("seq_batch", "sort_batch"))

FindClusterFreq(diagnosis, c("patient_id", "prognosis"), "sct_clusters") %>%
  group_by(cluster) %>%
  summarise(pval = wilcox.test(freq ~ prognosis)$p.value, .groups = "keep") %>%
  arrange(pval)

# Add Response Metadata
meta <- seurat[["patient_id"]] %>% rownames_to_column()

val <- read_excel(
  here("clinical_info/TARGET_AML_ClinicalData_Validation_20181213.xlsx")
)
dis <- read_excel(
  here("clinical_info/TARGET_AML_ClinicalData_Discovery_20181213.xlsx")
)
seq <- read_excel(file.path("../preprocessing/sample_info/Global Demultiplexing and Annotation.xlsx"))

clinical <- bind_rows(val, dis) %>%
  distinct() %>%
  separate(`TARGET USI`, into = c(NA, NA, "patient_id"), sep = "-")

patient_info <- filter(clinical, patient_id %in% pull(seq, patient_id))

seurat@misc[["patient_data"]] <- data

# patients we dont have clinical data for
# [1] "PAXMKT" "PAWAMB" "PARVEX" "PAUJMC" "PAUSFM" "PAUXYG" "PAXLWH" "PAWAMM" "PAVCNG" "PAUVMN" "PAUVEF" "PAUUWT"
# [13] "PAWZZN" "PAWUWJ" "PAUZSY" "PAVBFN" "PAWDTX" "PAUYAE" "PAUIIB" "PAUNVN"


## Frequency and Survival in scRNA seq Data ----
freq <- FindClusterFreq(diagnosis, c("patient_id", "prognosis"), "sct_clusters") %>%
  filter(cluster == 13) %>%
  dplyr::select(patient_id, freq)

summarise(freq,
  xbar = mean(freq),
  med = median(freq),
  q1 = quantile(freq)[2],
  q3 = quantile(freq)[4]
)

sc_survival_data <- patient_info %>%
  left_join(freq) %>%
  mutate(
    status = if_else(`First Event` == "Relapse", 1, 0),
    cluster_13_risk = if_else(freq >= 6.790474, "13 High", "13 Low")
  ) %>%
  dplyr::select(`Event Free Survival Time in Days`, status, cluster_13_risk, patient_id) %>%
  remove_missing()

sc_survival <- survfit(Surv(`Event Free Survival Time in Days`, status) ~ cluster_13_risk, data = sc_survival_data)
survdiff(Surv(`Event Free Survival Time in Days`, status) ~ cluster_13_risk, data = sc_survival_data)

ggsurvplot(
  data = sc_survival_data,
  fit = sc_survival,
  xlab = "Days",
  ylab = "Event Free Survival Probibility",
  risk.table = TRUE
)

## Create matrix for CIBERSORT with sc and TARGET Data ----

seurat_down <- subset(seurat, downsample = 100)

cibersort_data <- GetAssayData(seurat_down, slot = "data")
cibersort_data <- as.data.frame(cibersort_data)
cibersort_data2 <- rbind(paste0("cluster", Idents(seurat_down)), cibersort_data)

write.csv(cibersort_data2, file = here("cibersort.csv"), col.names = FALSE)

plan("multisession")

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ann <- getBM(
  attributes = c(
    "ensembl_gene_id_version",
    "hgnc_symbol"
  ),
  mart = ensembl
)

files <- list.files(here("target"), full.names = TRUE, recursive = TRUE) %>%
  map(data.table::fread)

ann <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "hgnc_symbol"
  ),
  mart = ensembl
)

plan("multisession")
res <- files %>%
  future_map(~ dplyr::select(., 1:2) %>%
    `colnames<-`(c("ensembl_gene_id", "data")) %>%
    filter(str_detect(ensembl_gene_id, pattern = "^ENS")) %>%
    right_join(ann) %>%
    as.tibble() %>%
    dplyr::select(hgnc_symbol, data) %>%
    group_by(hgnc_symbol) %>%
    summarise(data = sum(data))) %>%
  reduce(full_join, by = "hgnc_symbol") %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .)))

names <- list.files(here("target"), recursive = TRUE) %>%
  str_remove_all(".gene.quantification.txt")

colnames(res) <- c("GeneSymbol", names)

groups <- ceiling(ncol(res) / 500)
snames <- res[, 1]
res2 <- res[, -1]

for (b in 1:groups) {
  idx1 <- (b - 1) * 500
  idx2 <- b * 500
  if (idx2 > ncol(res2)) idx2 <- ncol(res2)
  res3 <- res2[-1, idx1:idx2]
  res3 <- cbind(snames[-1, ], res3)
  fname <- paste0("target_data_", b, ".txt")
  write_tsv(res3, file = here(fname))
}

## Import data from CIBERSORT and run survival ----
cibersort <- read_tsv(here("target_cibersort.txt")) %>%
  separate(col = "Mixture", into = c(NA, "numeric_id", "patient_id", "code1", "code2")) %>%
  inner_join(clinical)

summarise(cibersort,
  xbar = mean(cluster13),
  med = median(cluster13),
  q1 = quantile(cluster13)[2],
  q3 = quantile(cluster13)[4]
)

# survival with FLT3-ITD +
target_survival_data_itd <- cibersort %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "yes") %>%
  mutate(
    status = if_else(`First Event` == "relapse", 1, 0),
    cluster_13_risk = if_else(cluster13 >= mean(cluster13), "13 High", "13 Low")
  ) %>%
  dplyr::select(`Event Free Survival Time in Days`, status, cluster_13_risk, patient_id, `FLT3/ITD positive?`) %>%
  remove_missing()

target_ft3_survival <- survfit(Surv(`Event Free Survival Time in Days`, status) ~ cluster_13_risk, data = target_survival_data_itd)
survdiff(Surv(`Event Free Survival Time in Days`, status) ~ cluster_13_risk, data = target_survival_data_itd)

ggsurvplot(
  data = target_survival_data_itd,
  fit = target_ft3_survival,
  xlab = "Days",
  ylab = "Event Free Survival Probibility",
  risk.table = TRUE
)

# survival with FLT3-ITD-
target_survival_data_no_itd <- cibersort %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "no") %>%
  mutate(
    status = if_else(`First Event` == "relapse", 1, 0),
    cluster_13_risk = if_else(cluster13 >= mean(cluster13), "13 High", "13 Low")
  ) %>%
  dplyr::select(`Event Free Survival Time in Days`, status, cluster_13_risk, patient_id, `FLT3/ITD positive?`) %>%
  remove_missing()

target_no_ft3_survival <- survfit(Surv(`Event Free Survival Time in Days`, status) ~ cluster_13_risk, data = target_survival_data_no_itd)
survdiff(Surv(`Event Free Survival Time in Days`, status) ~ cluster_13_risk, data = target_survival_data_no_itd)

ggsurvplot(
  data = target_survival_data_no_itd,
  fit = target_no_ft3_survival,
  xlab = "Days",
  ylab = "Event Free Survival Probibility",
  risk.table = TRUE
)

## Lasso regression model
library(glmnet)
lasso_data <- cibersort %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "yes") %>%
  dplyr::select(
    cluster13, `WBC at Diagnosis`, `Bone marrow leukemic blast percentage (%)`,
    `Peripheral blasts (%)`, `Cytogenetic Complexity`,
    `FLT3/ITD allelic ratio`, `Event Free Survival Time in Days`
  ) %>%
  mutate(across(.fns = as.numeric)) %>%
  remove_missing()

lasso_model <- cv.glmnet(dplyr::select(lasso_data, -7) %>% as.matrix(), pull(lasso_data, 7))

lasso_model_final <- glmnet(dplyr::select(lasso_data, -7) %>% as.matrix(), pull(lasso_data, 7), lambda = lasso_model$lambda.min)

coef(lasso_model_final)
# (Intercept)                                 289.485315
# cluster13                                 -1127.728367
# WBC at Diagnosis                             -1.136043
# Bone marrow leukemic blast percentage (%)     2.790492
# Peripheral blasts (%)                         6.400624
# Cytogenetic Complexity                      -94.185055
# FLT3/ITD allelic ratio                        .

score <- cibersort %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "yes") %>%
  dplyr::select(
    cluster13, `WBC at Diagnosis`, `Bone marrow leukemic blast percentage (%)`,
    `Peripheral blasts (%)`, `Cytogenetic Complexity`, `First Event`,
    `FLT3/ITD allelic ratio`, `Event Free Survival Time in Days`
  ) %>%
  mutate(
    status = if_else(`First Event` == "relapse", 1, 0),
    `First Event` = NULL,
    across(.fns = as.numeric)
  ) %>%
  remove_missing() %>%
  set_names(c(
    "cluster13", "wbc", "leukemic_blasts", "peripheral_blasts",
    "cyto_complex", "allelic_ratio", "efs_days", "status"
  )) %>%
  mutate(
    score = -1127.73 * cluster13 - 1.14 * wbc + 2.79 * leukemic_blasts + 6.40 * peripheral_blasts - 94.19 * cyto_complex + 289.49,
    score_bin = if_else(score >= median(score), "High", "Low"),
    score_bin2 = cut(score, quantiles, labels = FALSE) %>% replace_na(1)
  )

score_survival <- survfit(Surv(efs_days, status) ~ score_bin2, data = score)
survdiff(Surv(efs_days, status) ~ score_bin2, data = score)

ggsurvplot(
  data = score,
  fit = score_survival,
  xlab = "Days",
  ylab = "Event Free Survival Probibility",
  risk.table = TRUE
)
