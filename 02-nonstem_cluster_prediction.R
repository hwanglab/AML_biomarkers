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

FindClusterFreq(diagnosis[[]], c("patient_id", "prognosis"), "clusters") %>%
  group_by(cluster) %>%
  summarise(pval = wilcox.test(freq ~ prognosis)$p.value, .groups = "keep") %>%
  arrange(pval)

FindClusterFreq(diagnosis[[]], c("patient_id", "prognosis"), "sct_clusters") %>%
  group_by(cluster) %>%
  dplyr::select(patient_id, freq, cluster)

summarise(freq,
          xbar = mean(freq),
          med = median(freq),
          q1 = quantile(freq)[2],
          q3 = quantile(freq)[4])

wilcox_clusters <- c("17", "3", "29", "5", "23", "18", "33", "26")

ReturnDifferences <- . %>%
  #mutate(freq = scale(freq, center = FALSE)) %>%
  group_by(status) %>%
  summarise(freq_mean = mean(freq)) %>%
  ungroup() %>%
  summarise(value = max(freq_mean) - min(freq_mean),
            value2 = status[which.max(freq_mean)]) %>%
  mutate(value3 = if_else(value2 == 1, value, -1 * value)) %>%
           pull(value3)

survival_models <- freq %>%
  right_join(clinical) %>%
  mutate(status = if_else(`First Event` == "Relapse", 1, 0),
                cluster_risk = if_else(freq >= mean(freq), "High", "Low")) %>%
  dplyr::select(`Event Free Survival Time in Days`, status, cluster_risk, patient_id, freq) %>%
  remove_missing() %>%
  filter(!near(freq, mean(freq))) %>%
  nest() %>%
  mutate(survival = map(data, ~ survdiff(Surv(`Event Free Survival Time in Days`, status) ~ cluster_risk, data = .x))) %>%
  mutate(chi_sq = map_dbl(survival, ~ .x[["chisq"]]),
         p_val = map_dbl(survival, CalculatePValues),
         log_p = -log(p_val),
         freq = map_dbl(data, ~ mean(.x[["freq"]])),
         wilcox = if_else(cluster %in% wilcox_clusters, "wilcox_sig", "not_wilcox_sig"),
         chi_sig = if_else(p_val <= 0.05, "chi_sig", "not_chi_sig"),
         sd = map_dbl(data, ~ sd(.x[["freq"]])),
         mean_diff = map_dbl(data, ReturnDifferences)) %>%
  arrange(p_val)

pdf(file = here("plots/sc_cluster_survival_analysis.pdf"), width = 3, height = 5)
ggplot(data = survival_models, mapping = aes(x = mean_diff, y = log_p)) +
  geom_hline(yintercept = -log(0.05)) +
  geom_point() +
  theme_classic() +
  ggrepel::geom_label_repel(data = filter(survival_models, p_val <= 0.05), mapping = aes(label = cluster)) +
  xlab("Distance Between Poor (+) and Favorable (-)") +
  ylab("-log(P value)")
graphics.off()

# export clusters for CIBERSORT
seurat_down <- subset(diagnosis, downsample = 100)

cibersort_data <- GetAssayData(seurat_down, slot = "data")
cibersort_data <- as.data.frame(cibersort_data)
cibersort_data2 <- rbind(paste0("cluster", Idents(seurat_down)), cibersort_data)

write.table(cibersort_data2, file = here("cibersort_diagnosis_nonstem.txt"), quote = FALSE, sep = "\t")

## LASSO Model with CIBERSORT ----
cibersort <- read_tsv(here("cibersort_results/cibersort_nonstem_target_results.txt")) %>%
  separate(col = "Mixture", into = c(NA, NA, "patient_id")) %>%
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
  mutate(chi_sq = map_dbl(survival, ~ .x[["chisq"]]),
         p_val = map_dbl(survival, CalculatePValues)) %>%
  arrange(p_val)

lasso_data <- cibersort %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "yes")

lasso_model <- DoLassoModel(lasso_data,
                            `Event Free Survival Time in Days`,
                            exclude = c(`P-value`:RMSE,
                                        matches("^Year|^Age|^Overall"),
                                        where(is_character)))

coef(lasso_model)


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
  mutate(score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
         score_bin = if_else(score >= median(score), "High", "Low"))

DoSurvialAnalysis(score, score$`Event Free Survival Time in Days`, score$status, score$score)

score_survival <- survfit(Surv(`Event Free Survival Time in Days`, status) ~ score_bin, data = score)
survdiff(Surv(`Event Free Survival Time in Days`, status) ~ score_bin, data = score)
coxph(Surv(`Event Free Survival Time in Days`, status) ~ score, data = score)

pdf(here("18-23-cibersort_score_target.pdf"))
ggsurvplot(data = score,
           fit = score_survival, 
           xlab = "Days", 
           ylab = "Event Free Survival Probibility",
           risk.table = TRUE)
graphics.off()


## LASSO Model with features (dont use) ----
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


target_data <- target %>%
  transposeDF(rowname_col = "patient") %>%
  dplyr::select(patient, any_of(genes)) %>%
  mutate(across(.cols = !any_of("patient"), .fns = as.numeric)) %>%
  separate(col = "patient", into = c(NA, NA, "patient_id", NA, NA)) %>%
  inner_join(clinical) %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "yes") %>%
  mutate(`FLT3/ITD allelic ratio` = as.numeric(`FLT3/ITD allelic ratio`))

lasso_data_full <- target_data %>%
  mutate(across(.fns = as.numeric)) %>%
  remove_missing()
  
lasso_model_genes <- DoLassoModel(lasso_data_full, `Event Free Survival Time in Days`)

genes_selected <- coef(lasso_model_genes) %>% as.data.frame() %>% filter(s0 != 0) %>% rownames()

lasso_data_reduced <- lasso_data_full %>%
  dplyr::select(any_of(genes_selected), `WBC at Diagnosis`, `Bone marrow leukemic blast percentage (%)`,
                           `Peripheral blasts (%)`, `Cytogenetic Complexity`,
                           `FLT3/ITD allelic ratio`, `Event Free Survival Time in Days`)
lasso_model_selected <- cv.glmnet(dplyr::select(lasso_data_reduced, -`Event Free Survival Time in Days`) %>% as.matrix, pull(lasso_data_reduced, `Event Free Survival Time in Days`))

lasso_model_selected_final <- glmnet(dplyr::select(lasso_data_reduced, -`Event Free Survival Time in Days`) %>% as.matrix, pull(lasso_data_reduced, `Event Free Survival Time in Days`), lambda = lasso_model_selected$lambda.min)
coef(lasso_model_selected_final)

score_features <- target_data %>%
  mutate(score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
         score_bin = if_else(score >= mean(score), "High", "Low"),
         status = if_else(`First Event` == "relapse", 1, 0)) %>%
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

## LASSO Model using Features ----
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

### Train on TARGET (training set) ----
target_model <- DoLassoModel(target_split$train,
             outcome = `Event Free Survival Time in Days`,
             exclude = `First Event`)

#test function
DoLassoModel(target_split$train,
             outcome = `Event Free Survival Time in Days`,
             exclude = c(`First Event`, AAR2, `FLT3/ITD allelic ratio`))

coef(target_model)

target_train_score <- target_split$train %>%
  mutate(score = UseLASSOModelCoefs(.data, coef(target_model)),
         score_bin = if_else(score >= median(score), "High", "Low"),
         status = if_else(`First Event` == "relapse", 1, 0)) %>%
  dplyr::select(score, score_bin, `Event Free Survival Time in Days`, status)

DoSurvialAnalysis(target_train_score,
                  time = `Event Free Survival Time in Days`,
                  status = status,
                  predictor = score,
                  group_by = score_bin,
                  description = "TARGET Training Data",
                  lasso = target_model)

## Test with TARGET (Test) ----
target_test_score <- target_split$test %>%
  mutate(score = UseLASSOModelCoefs(.data, coef(target_model)),
         score_bin = if_else(score >= median(score), "High", "Low"),
         status = if_else(`First Event` == "relapse", 1, 0)) %>%
  dplyr::select(score, score_bin, `Event Free Survival Time in Days`, status)

DoSurvialAnalysis(target_test_score,
                  time = `Event Free Survival Time in Days`,
                  status = status,
                  predictor = score,
                  group_by = score_bin,
                  description = "TARGET Test Data",
                  lasso = target_model)



### Test with TCGA on TARGET model ----
tcga_data <- read_tsv(here("tcga_cibersort.tsv"))
tcga_data <- transposeDF(tcga_data)

annotated_tcga <- read_tsv(here("tcga/clinical.cart.2021-04-22/clinical.tsv"), na = "'--") %>%
  right_join(tcga_data) %>%
  mutate(annotated_tcga, 
         score = UseLASSOModelCoefs(cur_data(), coef(target_model)),
         score_bin = if_else(score >= median(score), "High", "Low"),
         status = if_else(vital_status == "Dead", 1, 0))

DoSurvialAnalysis(annotated_tcga,
                  time = days_to_death,
                  status = status,
                  predictor = score,
                  group_by = score_bin,
                  description = "TCGA Data",
                  lasso = target_model)

### Test with BeatAML on TARGET model ----

beat_aml <- read_tsv(here("aml_ohsu_2018/data_RNA_Seq_expression_cpm.txt"))
beat_aml_clinical <- read_tsv(here("aml_ohsu_2018/data_clinical_sample.txt"), skip = 4)
beat_aml_survival <- read_tsv(here("aml_ohsu_2018/KM_Plot__Overall_Survival__(months).txt"))
beat_aml_clinical2 <- inner_join(beat_aml_clinical, beat_aml_survival, by = c("PATIENT_ID" = "Patient ID")) %>%
  mutate(status = if_else(SAMPLE_TIMEPOINT == "Relapse", 1, 0)) %>%
  dplyr::select(SAMPLE_ID, FLT3_ITD_CONSENSUS_CALL, OS_MONTHS, PB_BLAST_PERCENTAGE, status) %>%
  rename(`Peripheral blasts (%)` = PB_BLAST_PERCENTAGE)

beat_aml_data <- beat_aml %>%
  dplyr::select(-Entrez_Gene_Id) %>%
  transposeDF(rowname_col = SAMPLE_ID) %>%
  mutate(across(.cols = !starts_with("SAMPLE_ID"), .fns = as.numeric)) %>%
  inner_join(beat_aml_clinical2) %>%
  filter(FLT3_ITD_CONSENSUS_CALL == "Positive") %>%
  mutate(across(.cols = !any_of(c("FLT3_ITD_CONSENSUS_CALL", "SAMPLE_ID")), .fns = as.numeric),
         OS_DAYS = OS_MONTHS * 365.25 / 12,
         score = UseLASSOModelCoefs(cur_data(), coef(target_model)),
         score_bin = if_else(beat_aml_data$score >= mean(beat_aml_data$score), "High", "Low"))

DoSurvialAnalysis(beat_aml_data,
                  time = OS_DAYS,
                  status = status,
                  predictor = score,
                  group_by = score_bin,
                  description = "Beat-AML Data",
                  lasso = target_model)
