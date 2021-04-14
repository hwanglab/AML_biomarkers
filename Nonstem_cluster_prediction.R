library(Seurat)
library(SeuratDisk)
library(here)
library(rsinglecell)
library(readxl)
library(biomaRt)
library(survival)
library(survminer)
library(scorecard)
library(glmnet)
library(tidyverse)
library(furrr)

seurat <- LoadH5Seurat(file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat"),
                       assays = c("ADT", "SCT"))

meta <- seurat[["patient_id"]] %>% rownames_to_column()
val <- read_excel(here("clinical_info/TARGET_AML_ClinicalData_Validation_20181213.xlsx"))
dis <- read_excel(here("clinical_info/TARGET_AML_ClinicalData_Discovery_20181213.xlsx"))
seq <- read_excel(file.path("../preprocessing/sample_info/Global Demultiplexing and Annotation.xlsx"))

UseLASSOModelCoefs <- function(x, coef) {
  mod2 <- coef[-1, ] %>% as.data.frame()
  w_vec <- list(0)
  mod2 <- mod2[which(mod2[1] != 0), ,drop = FALSE]
  
  for (i in seq_along(rownames(mod3))) {
    col <- rownames(mod3)[i]
    vec <- x[[col]]
    w <- mod2[i, 1]
    w_vec[[i]] <- vec * w
  }
  res <- w_vec %>% purrr::reduce(rbind) %>% colSums() + mod[1, ]
  return(res)
}

clinical <- bind_rows(val, dis) %>%
  distinct() %>%
  separate(`TARGET USI`, into = c(NA, NA, "patient_id"), sep = "-")

data <- filter(clinical, patient_id %in% pull(seq, patient_id))

seurat@misc[["patient_data"]] <- data
seurat@misc[["target_data"]] <- clinical

diagnosis <- subset(seurat, timepoint == "Diagnosis" & stemness == "Nonstem")
diagnosis <- DoDimensionReductions(diagnosis, batch_vars = c("seq_batch", "sort_batch"))

diagnosis[["clusters"]] <- Idents(diagnosis)

saveRDS(diagnosis, file = here("seurat_diagnosis_nonstem.Rds"))

FindClusterFreq(diagnosis, c("patient_id", "prognosis"), "clusters") %>%
  group_by(cluster) %>%
  summarise(pval = wilcox.test(freq ~ prognosis)$p.value, .groups = "keep") %>%
  arrange(pval)

freq <- FindClusterFreq(diagnosis, c("patient_id", "prognosis"), "sct_clusters") %>%
  filter(cluster %in% c("17", "3", "29", "5", "23", "18", "33", "26")) %>%
  group_by(cluster) %>%
  dplyr::select(patient_id, freq, cluster)

summarise(freq,
          xbar = mean(freq),
          med = median(freq),
          q1 = quantile(freq)[2],
          q3 = quantile(freq)[4])

survival_models <- freq %>%
  right_join(patient_info) %>%
  mutate(status = if_else(`First Event` == "Relapse", 1, 0),
                cluster_risk = if_else(freq >= mean(freq), "High", "Low")) %>%
  dplyr::select(`Event Free Survival Time in Days`, status, cluster_risk, patient_id) %>%
  remove_missing() %>%
  nest() %>%
  mutate(survival = map(data, ~ survdiff(Surv(`Event Free Survival Time in Days`, status) ~ cluster_risk, data = .x))) %>%
  mutate(chi_sq = map_dbl(survival, ~ .x[["chisq"]])) %>%
  arrange(desc(chi_sq))

# export clusters for CIBERSORT
seurat_down <- subset(diagnosis, downsample = 100)

cibersort_data <- GetAssayData(seurat_down, slot = "data")
cibersort_data <- as.data.frame(cibersort_data)
cibersort_data2 <- rbind(paste0("cluster", Idents(seurat_down)), cibersort_data)

write.table(cibersort_data2, file = here("cibersort_diagnosis_nonstem.txt"), quote = FALSE, sep = "\t")

# import cibersort data
cibersort <- read_tsv(here("cibersort_results/cibersort_nonstem_target_results.txt")) %>%
  separate(col = "Mixture", into = c(NA, "numeric_id", "patient_id", "code1", "code2")) %>%
  inner_join(clinical)

summarise(cibersort, across(.cols = c("cluster18", "cluster23", "cluster5"),
                            .fns = c(x = mean, med = median)))

cibersort %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "yes") %>%
  mutate(status = if_else(`First Event` == "relapse", 1, 0),
         across(.cols = c("cluster18", "cluster23", "cluster5"),
                .fns = ~ if_else(.x >= mean(.x, na.rm = TRUE), "High", "Low"))) %>%
  dplyr::select(`Event Free Survival Time in Days`, status, cluster18,
                cluster23, cluster5, patient_id, `FLT3/ITD positive?`) %>%
  remove_missing() %>%
  pivot_longer(cols = starts_with("cluster")) %>%
  group_by(name) %>%
  nest() %>%
  mutate(survival = map(data, ~ survdiff(Surv(`Event Free Survival Time in Days`, status) ~ value, data = .x))) %>%
  mutate(chi_sq = map_dbl(survival, ~ .x[["chisq"]])) %>%
  arrange(desc(chi_sq))

lasso_data <- cibersort %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "yes") %>%
  dplyr::select(cluster18, cluster23, cluster5, `WBC at Diagnosis`,
                `Bone marrow leukemic blast percentage (%)`,
                `Peripheral blasts (%)`, `Cytogenetic Complexity`,
                `FLT3/ITD allelic ratio`, `Event Free Survival Time in Days`) %>%
  mutate(across(.fns = as.numeric)) %>%
  remove_missing()

lasso_model <- cv.glmnet(dplyr::select(lasso_data, -9) %>% as.matrix, pull(lasso_data, 9))

lasso_model_final <- glmnet(dplyr::select(lasso_data, -9) %>% as.matrix, pull(lasso_data, 9), lambda = lasso_model$lambda.min)

coef(lasso_model_final)

score <- cibersort %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "yes") %>%
  dplyr::select(cluster18, cluster23, cluster5, `WBC at Diagnosis`,
                `Bone marrow leukemic blast percentage (%)`,
                `Peripheral blasts (%)`, `Cytogenetic Complexity`, `First Event`,
                `FLT3/ITD allelic ratio`, `Event Free Survival Time in Days`) %>%
  mutate(status = if_else(`First Event` == "relapse", 1, 0),
         `First Event` = NULL,
         across(.fns = as.numeric)) %>%
  remove_missing() %>%
  set_names(c("cluster18", "cluster23", "cluster5", "wbc", "leukemic_blasts", 
              "peripheral_blasts", "cyto_complex", "allelic_ratio",
              "efs_days", "status")) %>%
  mutate(score = 1543.36 * cluster18 - 8063.56 * cluster23 - 2.15 * wbc + 3.44 * leukemic_blasts + 11.10 * peripheral_blasts - 118.07 * allelic_ratio + 256.39,
         score_bin = if_else(score >= median(score), "High", "Low"))

score_survival <- survfit(Surv(efs_days, status) ~ score_bin, data = score)
survdiff(Surv(efs_days, status) ~ score_bin, data = score)

pdf(here("18-23-cibersort_score_target.pdf"))
ggsurvplot(data = score,
           fit = score_survival, 
           xlab = "Days", 
           ylab = "Event Free Survival Probibility",
           risk.table = TRUE)
graphics.off()

gep <- read_tsv(here("cibersort_results/cibersort_nonstem_gep.txt"))

genes <- gep %>%
  dplyr::select(cluster18, cluster23, NAME) %>%
  filter(cluster18 != 1 | cluster23 != 1) %>%
  filter(cluster18 == 1 | cluster23 == 1) %>%
  pull(NAME)

target <- read_tsv(file = here("target_data.txt"))

clinical2 <- dplyr::select(clinical, `WBC at Diagnosis`, `Bone marrow leukemic blast percentage (%)`,
                    `Peripheral blasts (%)`, `Cytogenetic Complexity`,
                    `FLT3/ITD allelic ratio`, `Event Free Survival Time in Days`, `FLT3/ITD positive?`, patient_id)

lasso_data_full <- data.table::transpose(target) %>%
  slice(2:nrow(.)) %>%
  as_tibble() %>%
  set_names(target$GeneSymbol) %>%
  mutate(patient = colnames(target[2:ncol(target)]), .before = 1) %>%
  dplyr::select(patient, any_of(genes)) %>%
  mutate(across(.cols = 2:ncol(.), .fns = as.numeric)) %>%
  separate(col = "patient", into = c(NA, NA, "patient_id", NA, NA)) %>%
  inner_join(clinical2) %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "yes") %>%
  dplyr::select(-`FLT3/ITD positive?`, -patient_id) %>%
  mutate(across(.fns = as.numeric)) %>%
  remove_missing()
  
lasso_model_genes <- cv.glmnet(dplyr::select(lasso_data_full, -`Event Free Survival Time in Days`) %>% as.matrix, pull(lasso_data_full, `Event Free Survival Time in Days`))

lasso_model_genes_final <- glmnet(dplyr::select(lasso_data_full, -`Event Free Survival Time in Days`) %>% as.matrix, pull(lasso_data_full, `Event Free Survival Time in Days`), lambda = lasso_model_genes$lambda.min)

genes_selected <- coef(lasso_model_genes_final) %>% as.data.frame() %>% filter(s0 != 0) %>% rownames()

lasso_data_reduced <- lasso_data_full %>%
  dplyr::select(any_of(genes_selected), `WBC at Diagnosis`, `Bone marrow leukemic blast percentage (%)`,
                           `Peripheral blasts (%)`, `Cytogenetic Complexity`,
                           `FLT3/ITD allelic ratio`, `Event Free Survival Time in Days`, )
lasso_model_selected <- cv.glmnet(dplyr::select(lasso_data_reduced, -`Event Free Survival Time in Days`) %>% as.matrix, pull(lasso_data_reduced, `Event Free Survival Time in Days`))

lasso_model_selected_final <- glmnet(dplyr::select(lasso_data_reduced, -`Event Free Survival Time in Days`) %>% as.matrix, pull(lasso_data_reduced, `Event Free Survival Time in Days`), lambda = lasso_model_selected$lambda.min)
coef(lasso_model_selected_final)

score_features <- data.table::transpose(target) %>%
  slice(2:nrow(.)) %>%
  as_tibble() %>%
  set_names(target$GeneSymbol) %>%
  mutate(patient = colnames(target[2:ncol(target)]), .before = 1) %>%
  dplyr::select(patient, any_of(genes)) %>%
  mutate(across(.cols = 2:ncol(.), .fns = as.numeric)) %>%
  separate(col = "patient", into = c(NA, NA, "patient_id", NA, NA)) %>%
  inner_join(clinical) %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "yes") %>%
  mutate(`FLT3/ITD allelic ratio` = as.numeric(`FLT3/ITD allelic ratio`)) %>%
  group_by(across(where(is_character))) %>%
  mutate(score = UseLASSOModelCoefs(cur_data(), coef(lasso_model_selected_final)),
         # need to ungroup here
         score_bin = if_else(score >= mean(score), "High", "Low"),
         status = if_else(`First Event` == "relapse", 1, 0)) %>%
  ungroup() %>%
  dplyr::select(score, score_bin, `Event Free Survival Time in Days`, status)

score_survival <- survfit(Surv(`Event Free Survival Time in Days`, status) ~ score_bin, data = score_features)
survdiff(Surv(`Event Free Survival Time in Days`, status) ~ score_bin, data = score_features)

pdf(here("18-23-features_score_target.pdf"))
ggsurvplot(data = score_features,
           fit = score_survival, 
           xlab = "Days", 
           ylab = "Event Free Survival Probibility",
           risk.table = TRUE)
graphics.off()

## Test model on Beat AML ----

beat_aml <- read_tsv(here("aml_ohsu_2018/data_RNA_Seq_expression_cpm.txt"))
beat_aml_clinical <- read_tsv(here("aml_ohsu_2018/data_clinical_sample.txt"), skip = 4)
beat_aml_survival <- read_tsv(here("aml_ohsu_2018/KM_Plot__Overall_Survival__(months).txt"))
beat_aml_clinical2 <- inner_join(beat_aml_clinical, beat_aml_survival, by = c("PATIENT_ID" = "Patient ID")) %>%
  mutate(status = if_else(SAMPLE_TIMEPOINT == "Relapse", 1, 0)) %>%
  dplyr::select(SAMPLE_ID, FLT3_ITD_CONSENSUS_CALL, OS_MONTHS, PB_BLAST_PERCENTAGE, status) %>%
  rename(`Peripheral blasts (%)` = PB_BLAST_PERCENTAGE)

beat_aml_data <- beat_aml %>%
  dplyr::select(-Entrez_Gene_Id) %>%
  data.table::transpose() %>%
  slice(2:nrow(.)) %>%
  as_tibble() %>%
  set_names(beat_aml$Hugo_Symbol) %>%
  mutate(SAMPLE_ID = colnames(beat_aml[3:ncol(beat_aml)]), .before = 1) %>%
  mutate(across(.cols = !starts_with("SAMPLE_ID"), .fns = as.numeric))%>%
  inner_join(beat_aml_clinical2) %>%
  filter(FLT3_ITD_CONSENSUS_CALL == "Positive") %>%
  dplyr::select(-FLT3_ITD_CONSENSUS_CALL, -SAMPLE_ID) %>%
  remove_missing() %>%
  mutate(across(.fns = as.numeric),
         OS_DAYS = OS_MONTHS * 365.25 / 12,
         score = UseLASSOModelCoefs(.data, coef(lasso_model_selected_final)),
         score_bin = if_else(score >= mean(score), "High", "Low"))


score_survival <- survfit(Surv(OS_DAYS, status) ~ score_bin, data = beat_aml_data)
survdiff(Surv(OS_DAYS, status) ~ score_bin, data = beat_aml_data)

pdf(here("plots/18-23-features_score_BeatAML.pdf"))
ggsurvplot(data = beat_aml_data,
           fit = score_survival, 
           xlab = "Days", 
           ylab = "Event Free Survival Probibility",
           risk.table = TRUE)
graphics.off()

# Test TARGET Without NOP53 and PCAT18
coef_for_beat_aml <- coef(lasso_model_selected_final)
coef_for_beat_aml <- coef_for_beat_aml[!(rownames(coef_for_beat_aml) %in% c("NOP53", "PCAT18")), , drop = FALSE]
score_features2 <- data.table::transpose(target) %>%
  slice(2:nrow(.)) %>%
  as_tibble() %>%
  set_names(target$GeneSymbol) %>%
  mutate(patient = colnames(target[2:ncol(target)]), .before = 1) %>%
  dplyr::select(patient, any_of(genes)) %>%
  mutate(across(.cols = 2:ncol(.), .fns = as.numeric)) %>%
  separate(col = "patient", into = c(NA, NA, "patient_id", NA, NA)) %>%
  inner_join(clinical) %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "yes") %>%
  mutate(score = UseLASSOModelCoefs(.data, coef_for_beat_aml),
         score_bin = if_else(score >= median(score), "High", "Low"),
         status = if_else(`First Event` == "relapse", 1, 0)) %>%
  dplyr::select(score, score_bin, `Event Free Survival Time in Days`, status)

score_survival <- survfit(Surv(`Event Free Survival Time in Days`, status) ~ score_bin, data = score_features2)
survdiff(Surv(`Event Free Survival Time in Days`, status) ~ score_bin, data = score_features2)

## Split TARGET into training + test ----

clinical2 <- dplyr::select(clinical, `WBC at Diagnosis`, `Bone marrow leukemic blast percentage (%)`,
                           `Peripheral blasts (%)`, `Cytogenetic Complexity`,
                           `FLT3/ITD allelic ratio`, `Event Free Survival Time in Days`, `FLT3/ITD positive?`, patient_id, `First Event`)

target_before_split <- data.table::transpose(target) %>%
  slice(2:nrow(.)) %>%
  as_tibble() %>%
  set_names(target$GeneSymbol) %>%
  mutate(patient = colnames(target[2:ncol(target)]), .before = 1) %>%
  dplyr::select(patient, any_of(genes)) %>%
  mutate(across(.cols = 2:ncol(.), .fns = as.numeric)) %>%
  separate(col = "patient", into = c(NA, NA, "patient_id", NA, NA)) %>%
  inner_join(clinical2) %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "yes") %>%
  dplyr::select(-`FLT3/ITD positive?`, -patient_id) %>%
  mutate(across(.cols = !`First Event`, .fns = as.numeric)) %>%
  remove_missing()

target_split <- split_df(target_before_split, y = "Event Free Survival Time in Days", ratios = c(0.7, 0.3))

lasso_model_genes <- cv.glmnet(dplyr::select(target_split$train, -`Event Free Survival Time in Days`, -`First Event`) %>% as.matrix, pull(target_split$train, `Event Free Survival Time in Days`))

lasso_model_genes_final <- glmnet(dplyr::select(target_split$train, -`Event Free Survival Time in Days`, -`First Event`) %>% as.matrix, pull(target_split$train, `Event Free Survival Time in Days`), lambda = lasso_model_genes$lambda.min)

genes_selected <- coef(lasso_model_genes_final) %>% as.data.frame() %>% filter(s0 != 0) %>% rownames()

lasso_data_reduced <- target_split$train %>%
  dplyr::select(any_of(genes_selected), `WBC at Diagnosis`, `Bone marrow leukemic blast percentage (%)`,
                `Peripheral blasts (%)`, `Cytogenetic Complexity`,
                `FLT3/ITD allelic ratio`, `Event Free Survival Time in Days`)
lasso_model_selected <- cv.glmnet(dplyr::select(lasso_data_reduced, -`Event Free Survival Time in Days`) %>% as.matrix, pull(lasso_data_reduced, `Event Free Survival Time in Days`))

lasso_model_selected_final <- glmnet(dplyr::select(lasso_data_reduced, -`Event Free Survival Time in Days`) %>% as.matrix, pull(lasso_data_reduced, `Event Free Survival Time in Days`), lambda = lasso_model_selected$lambda.min)
coef(lasso_model_selected_final)

# test on test cohort
target_test_score <- target_split$test %>%
  mutate(score = UseLASSOModelCoefs(.data, coef(lasso_model_selected_final)),
         score_bin = if_else(score >= median(score), "High", "Low"),
         status = if_else(`First Event` == "relapse", 1, 0)) %>%
  dplyr::select(score, score_bin, `Event Free Survival Time in Days`, status)

score_survival <- survfit(Surv(`Event Free Survival Time in Days`, status) ~ score_bin, data = target_test_score)
survdiff(Surv(`Event Free Survival Time in Days`, status) ~ score_bin, data = target_test_score)

ggsurvplot(data = target_test_score,
           fit = score_survival, 
           xlab = "Days", 
           ylab = "Event Free Survival Probibility",
           risk.table = TRUE)

# plots for training cohort
target_train_score <- target_split$train %>%
  mutate(score = UseLASSOModelCoefs(.data, coef(lasso_model_selected_final)),
         score_bin = if_else(score >= median(score), "High", "Low"),
         status = if_else(`First Event` == "relapse", 1, 0)) %>%
  dplyr::select(score, score_bin, `Event Free Survival Time in Days`, status)

score_survival_train <- survfit(Surv(`Event Free Survival Time in Days`, status) ~ score_bin, data = target_train_score)
survdiff(Surv(`Event Free Survival Time in Days`, status) ~ score_bin, data = target_train_score)

ggsurvplot(data = target_train_score,
           fit = score_survival_train, 
           xlab = "Days", 
           ylab = "Event Free Survival Probibility",
           risk.table = TRUE)
