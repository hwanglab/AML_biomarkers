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
library(xfun)
library(msigdbr)
library(furrr)
library(EnhancedVolcano)
library(tidyverse)

source(here("functions.R"))

val <- read_excel(here("clinical_info/TARGET_AML_ClinicalData_Validation_20181213.xlsx"))
dis <- read_excel(here("clinical_info/TARGET_AML_ClinicalData_Discovery_20181213.xlsx"))
seq <- read_excel(file.path("../preprocessing/sample_info/Global Demultiplexing and Annotation.xlsx"))

clinical <- bind_rows(val, dis) %>%
  distinct() %>%
  separate(`TARGET USI`, into = c(NA, NA, "patient_id"), sep = "-")

sig_level <- 0.10
rerun_cibersort <- FALSE
rerun_cellphone <- FALSE

## Subset + Dimension Reductions ----
diagnosis <- cache_rds(expr = {
  seurat <- LoadH5Seurat(file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat"),
                        assays = c("ADT", "SCT"))
  
  DefaultAssay(seurat) <- "SCT"

  meta <- seurat[["patient_id"]] %>% rownames_to_column()

  clinical <- bind_rows(val, dis) %>%
    distinct() %>%
    separate(`TARGET USI`, into = c(NA, NA, "patient_id"), sep = "-")

  data <- filter(clinical, patient_id %in% pull(seq, patient_id))

  seurat@misc[["patient_data"]] <- data
  seurat@misc[["target_data"]] <- clinical

  diagnosis <- subset(seurat, timepoint == "Diagnosis" & stemness == "Nonstem")
  diagnosis <- DoDimensionReductions(diagnosis, batch_vars = c("seq_batch", "sort_batch"))

  diagnosis[["clusters"]] <- Idents(diagnosis)
  diagnosis
  },
  file = "02-seurat_diagnosis_nonstem.rds")

## Find DE Clusters ----
wilcox_clusters <- FindClusterFreq(diagnosis[[]], c("patient_id", "prognosis"), "clusters") %>%
  group_by(cluster) %>%
  summarise(pval = wilcox.test(freq ~ prognosis)$p.value, .groups = "keep") %>%
  filter(pval <= sig_level) %>%
  pull(cluster)

freq <- FindClusterFreq(diagnosis[[]], c("patient_id", "prognosis"), "clusters") %>%
  group_by(cluster) %>%
  dplyr::select(patient_id, freq, cluster)

summarise(freq,
          xbar = mean(freq),
          med = median(freq),
          q1 = quantile(freq)[2],
          q3 = quantile(freq)[4])

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
         p_val = map_dbl(survival, CalculatePValues.chi),
         log_p = -log(p_val),
         freq = map_dbl(data, ~ mean(.x[["freq"]])),
         wilcox = if_else(cluster %in% wilcox_clusters, "wilcox_sig", "not_wilcox_sig"),
         chi_sig = if_else(p_val <= sig_level, "chi_sig", "not_chi_sig"),
         sd = map_dbl(data, ~ sd(.x[["freq"]])),
         mean_diff = map_dbl(data, ReturnDifferences)) %>%
  arrange(p_val)

pdf(file = here("plots/sc_cluster_survival_analysis.pdf"), width = 3, height = 5)
ggplot(data = survival_models, mapping = aes(x = mean_diff, y = log_p)) +
  geom_hline(yintercept = -log(sig_level)) +
  geom_point() +
  theme_classic() +
  ggrepel::geom_label_repel(data = filter(survival_models, p_val <= sig_level), mapping = aes(label = cluster)) +
  xlab("Distance Between Poor (+) and Favorable (-)") +
  ylab("-log(P value)")
graphics.off()

## Run GSVA on Clusters ----

gene_sets <- list(HALLMARK = msigdbr(species = "Homo sapiens", category = "H"),
                  GO       = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP"),
                  REACTOME = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME"),
                  KEGG     = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG"),
                  ONCO     = msigdbr(species = "Homo sapiens", category = "C6")) %>%
  map(~ ListGeneSets(.x, gene_set = "name", gene_name = "symbol"))

gsva_res <- map(gene_sets, ~ RunGSVA(diagnosis, gene_sets = .x, replicates = 3))

stat_res <- map2(gsva_res, names(gsva_res), RunStats, p = 0.05)

pdf(file = here("plots/HALLMARK_GSVA_volcano_plots.pdf"), onefile = TRUE)
stat_res$HALLMARK$plots
graphics.off()

pdf(file = here("plots/GO_GSVA_volcano_plots.pdf"), onefile = TRUE)
stat_res$GO$plots
graphics.off()

pdf(file = here("plots/REACTOME_GSVA_volcano_plots.pdf"), onefile = TRUE)
stat_res$REACTOME$plots
graphics.off()

pdf(file = here("plots/KEGG_GSVA_volcano_plots.pdf"), onefile = TRUE)
stat_res$KEGG$plots
graphics.off()

pdf(file = here("plots/ONCO_GSVA_volcano_plots.pdf"), onefile = TRUE)
stat_res$ONCO$plots
graphics.off()

top_tables <- map(stat_res, ~ .x[["Top Table"]])

top_tables %>%
  reduce(bind_rows) %>%
  list("All Gene Sets" = .) %>%
  append(top_tables) %>%
  openxlsx::write.xlsx(file = here("outs/GSVA_DE_results.xlsx"))

## Run CellPhoneDB ----

if (rerun_cellphone) {
  seurat_tmp <- LoadH5Seurat(file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat"),
                         assays = c("RNA"))
  DefaultAssay(seurat_tmp) <- "RNA"
  tmp <- subset(seurat_tmp, timepoint == "Diagnosis" & stemness == "Nonstem")
  tmp <- NormalizeData(tmp)
  
  Idents(diagnosis) %>%
    as_tibble(rownames = "Cell") %>%
    rename(cell_type = value) %>%
    write_tsv(here("data/cellphonedb_in/metadata.tsv"))
  
  GetAssayData(tmp, "data") %>%
    DropletUtils::write10xCounts(path = here("data/cellphonedb_in/data"))
  
  file.rename(here("data/cellphonedb_in/data/genes.tsv"),
              here("data/cellphonedb_in/data/features.tsv"))
  
  rm(tmp, seurat_tmp)
  
  system(here("cellphonedb.sh"))
}

## DE on Clusters ----
sig_clusters <- survival_models %>%
  filter(chi_sig == "chi_sig") %>%
  pull(cluster)

de_results <- map(sig_clusters, ~ FindMarkers(diagnosis,
                                              ident.1 = .x,
                                              test.use = "MAST"))
names(de_results) <- sig_clusters

de_name <- paste0(sig_clusters, collapse = "+")
de_results[[de_name]] <- FindMarkers(diagnosis,
                                     ident.1 = sig_clusters,
                                     test.use = "MAST")

de_results[["All Markers"]] <- FindAllMarkers(diagnosis, test.use = "MAST")
openxlsx::write.xlsx(de_results, 
                     file = here("outs/cluster_DE_results.xlsx"), 
                     rowNames = TRUE)

EnhancedVolcano(de_results$`20`, x = "avg_log2FC", y = "adj_p_val")

## Run CIBERSORTx ----
if (rerun_cibersort) {
  seurat_down <- subset(diagnosis, downsample = 100)

  cibersort_data <- GetAssayData(seurat_down, slot = "data")
  cibersort_data <- as.data.frame(cibersort_data)
  cibersort_data2 <- rbind(paste0("cluster", Idents(seurat_down)), cibersort_data)

  write.table(cibersort_data2, file = here("cibersort_in/nonstem_clusters.txt"), quote = FALSE, sep = "\t")

  system(here("CIBERSORTx.sh"))
}

# include non-flt3 pats

## LASSO Model with CIBERSORT ----
sig_clusters <- survival_models %>%
  filter(chi_sig == "chi_sig") %>%
  pull(cluster) %>%
  paste0("cluster", .)

target_ciber <- read_tsv(here("outs/cibersort_results/CIBERSORTx_target_Results.txt")) %>%
  separate(col = "Mixture", into = c(NA, NA, "patient_id")) %>%
  select(all_of(sig_clusters), patient_id) %>%
  inner_join(clinical)

summarise(target_ciber, across(.cols = sig_clusters,
                               .fns = c(x = mean, med = median)))

target_ciber %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "yes") %>%
  mutate(status = if_else(`First Event` == "relapse", 1, 0),
         across(.cols = starts_with("cluster"),
                .fns = ~ if_else(.x >= mean(.x, na.rm = TRUE), "High", "Low"))) %>%
  dplyr::select(`Event Free Survival Time in Days`, status, 
                all_of(sig_clusters), patient_id, `FLT3/ITD positive?`) %>%
  remove_missing() %>%
  pivot_longer(cols = starts_with("cluster")) %>%
  group_by(name) %>%
  nest() %>%
  mutate(survival = map(data, ~ survdiff(Surv(`Event Free Survival Time in Days`, status) ~ value, data = .x))) %>%
  mutate(chi_sq = map_dbl(survival, ~ .x[["chisq"]]),
         p_val = map_dbl(survival, CalculatePValues.chi)) %>%
  arrange(p_val) %>%
  select(name, chi_sq, p_val) %>%
  write_csv(file = here("outs/target_single_cluster_survival.csv"))

### Train on TARGET (train) ----

target_ciber_split <- target_ciber %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "yes") %>%
  split_df(y = "Event Free Survival Time in Days", ratios = c(0.7, 0.3))

lasso_model <- DoLassoModel(target_ciber_split$train,
                            `Event Free Survival Time in Days`,
                            exclude = c(matches("^Year|^Age|^Overall"),
                                        where(is_character)))

coef(lasso_model) %>%
  as.data.frame() %>%
  rownames_to_column(var = "predictor") %>%
  rename(coef = s0) %>%
  write_csv(here("outs/TARGET_CIBERSORTx_LASSO_model.csv"))

target_ciber_split$train <- target_ciber_split$train %>%
  as_tibble() %>%
  mutate(status = if_else(`First Event` == "relapse", 1, 0),
         `First Event` = NULL) %>%
  mutate(score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
         score_bin = if_else(score >= median(score), "High", "Low"))

target_train <- DoSurvialAnalysis(target_ciber_split$train,
                  `Event Free Survival Time in Days`,
                  status,
                  score,
                  group_by = score_bin,
                  description = "TARGET Training",
                  lasso = lasso_model)

target_train_os <- DoSurvialAnalysis(target_ciber_split$train,
                                     `Overall Survival Time in Days`,
                                     status,
                                     score,
                                     group_by = score_bin,
                                     description = "TARGET Training",
                                     lasso = lasso_model)

### Test on TARGET (test) ----
target_ciber_split$test <- target_ciber_split$test %>%
  as_tibble() %>%
  mutate(status = if_else(`First Event` == "relapse", 1, 0),
         `First Event` = NULL) %>%
  mutate(score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
         score_bin = if_else(score >= median(score), "High", "Low"))

target_test <- DoSurvialAnalysis(target_ciber_split$test,
                  `Event Free Survival Time in Days`,
                  status,
                  score,
                  group_by = score_bin,
                  description = "TARGET Test",
                  lasso = lasso_model)

target_test_os <- DoSurvialAnalysis(target_ciber_split$test,
                                    `Overall Survival Time in Days`,
                                    status,
                                    score,
                                    group_by = score_bin,
                                    description = "TARGET Training",
                                    lasso = lasso_model)

### Test on TARGET (FLT3-ITD Negative) ----
target_ciber_flt_neg <- target_ciber %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "no")

target_ciber_flt3_neg <- target_ciber_flt_neg %>%
  mutate(status = if_else(`First Event` == "relapse", 1, 0),
         `First Event` = NULL) %>%
  mutate(score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
         score_bin = if_else(score >= median(score), "High", "Low"))

target_flt3 <- DoSurvialAnalysis(target_ciber_flt3_neg,
                                 `Event Free Survival Time in Days`,
                                 status,
                                 score,
                                 group_by = score_bin,
                                 description = "TARGET Test",
                                 lasso = lasso_model)

target_flt3_os <- DoSurvialAnalysis(target_ciber_flt3_neg,
                                    `Overall Survival Time in Days`,
                                    status,
                                    score,
                                    group_by = score_bin,
                                    description = "TARGET Training",
                                    lasso = lasso_model)


### Test on TCGA ----

tcga_ciber <- read_tsv(here("outs/cibersort_results/CIBERSORTx_tcga_Results.txt"))
tcga_ann <- read_tsv(here("data/tcga/clinical.cart.2021-04-22/clinical.tsv"), na = "'--")

tcga_ciber <- tcga_ciber %>%
  left_join(tcga_ann2, by = c("Mixture" = "case_submitter_id")) %>%
  filter(flt3_status == "Positive") %>%
  mutate(score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
         score_bin = if_else(score >= median(score), "High", "Low"),
         status = if_else(vital_status == "Dead", 1, 0),
         days_to_death = as.numeric(days_to_death))

tcga_surv <- DoSurvialAnalysis(tcga_ciber,
                                 days_to_death,
                                 status,
                                 score,
                                 group_by = score_bin,
                                 description = "TCGA (50% features present)",
                                 lasso = lasso_model)

### Test on BeatAML ----

beat_aml_clinical <- read_tsv(here("data/aml_ohsu_2018/data_clinical_sample.txt"), skip = 4)
beat_aml_survival <- read_tsv(here("data/aml_ohsu_2018/KM_Plot__Overall_Survival__(months).txt"))
beat_aml_clinical2 <- inner_join(beat_aml_clinical, beat_aml_survival, by = c("PATIENT_ID" = "Patient ID")) %>%
  mutate(status = if_else(SAMPLE_TIMEPOINT == "Relapse", 1, 0)) %>%
  dplyr::select(SAMPLE_ID, FLT3_ITD_CONSENSUS_CALL, OS_MONTHS, PB_BLAST_PERCENTAGE, status) %>%
  rename(`Peripheral blasts (%)` = PB_BLAST_PERCENTAGE)

beat_aml_ciber <- read_tsv(here("outs/cibersort_results/CIBERSORTx_beatAML_Results.txt")) %>%
  left_join(beat_aml_clinical2, by = c("Mixture" = "SAMPLE_ID")) %>%
  filter(FLT3_ITD_CONSENSUS_CALL == "Positive") %>%
  mutate(score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
         score_bin = if_else(score >= mean(score, na.rm = TRUE), "High", "Low"),
         OS_days = OS_MONTHS * (365.25 / 12))

beat_surv <- DoSurvialAnalysis(beat_aml_ciber,
                               OS_days,
                               status,
                               score,
                               group_by = score_bin,
                               description = "BeatAML (80% features present)",
                               lasso = lasso_model)

### Print Plots ----

pdf(file = here("plots/CIBERSORT_survival.pdf"))
target_train$plot
target_test$plot
target_train_os$plot
target_test_os$plot
target_flt3$plot
target_flt3_os$plot
tcga_surv$plot
beat_surv$plot
graphics.off()

## LASSO Model with ONLY CIBERSORTx features ----

target_ciber_no_filter <- read_tsv(here("outs/cibersort_results/CIBERSORTx_target_Results.txt")) %>%
  separate(col = "Mixture", into = c(NA, NA, "patient_id")) %>%
  inner_join(clinical)

target_ciber_no_filter %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "yes") %>%
  mutate(status = if_else(`First Event` == "relapse", 1, 0),
         across(.cols = starts_with("cluster"),
                .fns = ~ if_else(.x >= mean(.x, na.rm = TRUE), "High", "Low"))) %>%
  dplyr::select(`Event Free Survival Time in Days`, status, 
                starts_with("cluster"), patient_id, `FLT3/ITD positive?`) %>%
  remove_missing() %>%
  pivot_longer(cols = starts_with("cluster")) %>%
  group_by(name) %>%
  nest() %>%
  mutate(survival = map(data, ~ survdiff(Surv(`Event Free Survival Time in Days`, status) ~ value, data = .x))) %>%
  mutate(chi_sq = map_dbl(survival, ~ .x[["chisq"]]),
         p_val = map_dbl(survival, CalculatePValues.chi)) %>%
  arrange(p_val) %>%
  select(name, chi_sq, p_val) %>%
  write_csv(file = here("outs/target_single_cluster_survival_CIBERSORTx_only.csv"))

### Train on TARGET (train) ----

target_ciber_nf_split <- target_ciber_no_filter %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "yes") %>%
  split_df(y = "Event Free Survival Time in Days", ratios = c(0.7, 0.3))

lasso_model_nf <- DoLassoModel(target_ciber_nf_split$train,
                            `Event Free Survival Time in Days`,
                            exclude = c(matches("^Year|^Age|^Overall|^WBC|%"),
                                        where(is_character),
                                        `P-value`:RMSE))

coef(lasso_model_nf) %>%
  as.data.frame() %>%
  rownames_to_column(var = "predictor") %>%
  rename(coef = s0) %>%
  write_csv(here("outs/TARGET_CIBERSORTx_LASSO_model_only_CIBERSORTx.csv"))

target_ciber_nf_split$train <- target_ciber_nf_split$train %>%
  as_tibble() %>%
  mutate(status = if_else(`First Event` == "relapse", 1, 0),
         `First Event` = NULL) %>%
  mutate(score = UseLASSOModelCoefs(cur_data(), coef(lasso_model_nf)),
         score_bin = if_else(score >= median(score), "High", "Low"))

target_train_nf <- DoSurvialAnalysis(target_ciber_split$train,
                                  `Event Free Survival Time in Days`,
                                  status,
                                  score,
                                  group_by = score_bin,
                                  description = "TARGET Training",
                                  lasso = lasso_model_nf)

target_train_os_nf <- DoSurvialAnalysis(target_ciber_nf_split$train,
                                        `Overall Survival Time in Days`,
                                        status,
                                        score,
                                        group_by = score_bin,
                                        description = "TARGET Training",
                                        lasso = lasso_model_nf)

### Test on TARGET (test) ----

target_ciber_nf_split$test <- target_ciber_nf_split$test %>%
  as_tibble() %>%
  mutate(status = if_else(`First Event` == "relapse", 1, 0),
         `First Event` = NULL) %>%
  mutate(score = UseLASSOModelCoefs(cur_data(), coef(lasso_model_nf)),
         score_bin = if_else(score >= median(score), "High", "Low"))

target_test_nf <- DoSurvialAnalysis(target_ciber_nf_split$test,
                                 `Event Free Survival Time in Days`,
                                 status,
                                 score,
                                 group_by = score_bin,
                                 description = "TARGET Test",
                                 lasso = lasso_model_nf)

target_test_os_nf <- DoSurvialAnalysis(target_ciber_nf_split$test,
                                    `Overall Survival Time in Days`,
                                    status,
                                    score,
                                    group_by = score_bin,
                                    description = "TARGET Training",
                                    lasso = lasso_model_nf)

### Test on TCGA ----

tcga_ciber <- read_tsv(here("outs/cibersort_results/CIBERSORTx_tcga_Results.txt"))
tcga_ann <- read_tsv(here("data/tcga/clinical.cart.2021-04-22/clinical.tsv"), na = "'--")

tcga_ciber <- tcga_ciber %>%
  left_join(tcga_ann2, by = c("Mixture" = "case_submitter_id")) %>%
  filter(flt3_status == "Positive") %>%
  mutate(score = UseLASSOModelCoefs(cur_data(), coef(lasso_model_nf)),
         score_bin = if_else(score >= median(score), "High", "Low"),
         status = if_else(vital_status == "Dead", 1, 0),
         days_to_death = as.numeric(days_to_death))

tcga_surv_nf <- DoSurvialAnalysis(tcga_ciber,
                               days_to_death,
                               status,
                               score,
                               group_by = score_bin,
                               description = "TCGA (33% features present)",
                               lasso = lasso_model_nf)

### Test on BeatAML ----
beat_aml_ciber_nf <- read_tsv(here("outs/cibersort_results/CIBERSORTx_beatAML_Results.txt")) %>%
  left_join(beat_aml_clinical2, by = c("Mixture" = "SAMPLE_ID")) %>%
  filter(FLT3_ITD_CONSENSUS_CALL == "Positive") %>%
  mutate(score = UseLASSOModelCoefs(cur_data(), coef(lasso_model_nf)),
         score_bin = if_else(score >= mean(score, na.rm = TRUE), "High", "Low"),
         OS_days = OS_MONTHS * (365.25 / 12))

beat_surv_nf <- DoSurvialAnalysis(beat_aml_ciber_nf,
                               OS_days,
                               status,
                               score,
                               group_by = score_bin,
                               description = "BeatAML (67% features present)",
                               lasso = lasso_model)

### Print Plots ----
pdf(file = here("plots/CIBERSORT_survival_no_filtering.pdf"))
target_train_nf$plot
target_test_nf$plot
target_train_os_nf$plot
target_test_os_nf$plot
tcga_surv_nf$plot
beat_surv_nf$plot
graphics.off()

## LASSO Model with ONLY CIBERSORTx features ----

### Train on TCGA ----

tcga_ciber <- read_tsv(here("outs/cibersort_results/CIBERSORTx_tcga_Results.txt"))
tcga_ciber <- left_join(tcga_ciber, tcga_ann2, by = c("Mixture" = "case_submitter_id"))

lasso_model_nf <- DoLassoModel(tcga_ciber,
                               days_to_death,
                               exclude = c(matches("^Year|^Age|^Overall|^WBC|%|^lab_|^year_|^age_|cytogen|follow|birth|initial|diagnosis|result"),
                                           where(is_character),
                                           where(is_logical),
                                           where(lubridate::is.instant),
                                           `P-value`:RMSE,
                                           Mixture,
                                           patient_id,
                                           hydroxyurea_agent_administered_day_count,
                                           cumulative_agent_total_dose))

tcga_res <- tcga_ciber %>%
  rename(case_submitter_id = Mixture) %>%
  filter(flt3_status == "Positive") %>%
  mutate(score = UseLASSOModelCoefs(cur_data(), coef(lasso_model_nf)),
         score_bin = if_else(score >= median(score), "High", "Low"),
         status = if_else(vital_status == "Dead", 1, 0),
         days_to_death = as.numeric(days_to_death))

tcga_surv_nf <- DoSurvialAnalysis(tcga_res,
                                  days_to_death,
                                  status,
                                  score,
                                  group_by = score_bin,
                                  description = "TCGA",
                                  lasso = lasso_model_nf)

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
###   LASSO MODEL SIGNIFICANT   ###
###       FINISH ANALYSIS       ###
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

## LASSO Model using Features ----
clinical2 <- dplyr::select(clinical, `WBC at Diagnosis`, `Bone marrow leukemic blast percentage (%)`,
                           `Peripheral blasts (%)`, `Cytogenetic Complexity`,
                           `FLT3/ITD allelic ratio`, `Event Free Survival Time in Days`, `FLT3/ITD positive?`, patient_id, `First Event`,
                           `Overall Survival Time in Days`)

gep <- read_tsv(here("cibersort_in/nonstem_GEP.txt"))

genes <- gep %>%
  dplyr::select(cluster18, cluster23, genesymbols) %>%
  filter(cluster18 != 1 | cluster23 != 1) %>%
  filter(cluster18 == 1 | cluster23 == 1) %>%
  pull(genesymbols)

target <- read_tsv(file = here("cibersort_in/target_data.txt"))

target_before_split <- target %>%
  transposeDF() %>%
  dplyr::select(samples, any_of(genes)) %>%
  mutate(across(.cols = 2:ncol(.), .fns = as.numeric)) %>%
  separate(col = "samples", into = c(NA, NA, "patient_id", NA, NA)) %>%
  inner_join(clinical2) %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "yes") %>%
  dplyr::select(-`FLT3/ITD positive?`, -patient_id) %>%
  mutate(across(.cols = !`First Event`, .fns = as.numeric)) %>%
  remove_missing()

target_split <- split_df(target_before_split, y = "Event Free Survival Time in Days", ratios = c(0.7, 0.3))

### Train on TARGET (train) ----
target_model <- DoLassoModel(target_split$train,
             outcome = `Event Free Survival Time in Days`,
             exclude = c(`First Event`, `Overall Survival Time in Days`))

coef(target_model)

target_train_score <- target_split$train %>%
  mutate(score = UseLASSOModelCoefs(cur_data(), coef(target_model)),
         score_bin = if_else(score >= median(score), "High", "Low"),
         status = if_else(`First Event` == "relapse", 1, 0))

target_train_surv <- DoSurvialAnalysis(target_train_score,
                  time = `Event Free Survival Time in Days`,
                  status = status,
                  predictor = score,
                  group_by = score_bin,
                  description = "TARGET Training Data",
                  lasso = target_model)

target_train_surv_os <- DoSurvialAnalysis(target_train_score,
                                       time = `Overall Survival Time in Days`,
                                       status = status,
                                       predictor = score,
                                       group_by = score_bin,
                                       description = "TARGET Training Data",
                                       lasso = target_model)

### Test with TARGET (test) ----
target_test_score <- target_split$test %>%
  mutate(score = UseLASSOModelCoefs(cur_data(), coef(target_model)),
         score_bin = if_else(score >= median(score), "High", "Low"),
         status = if_else(`First Event` == "relapse", 1, 0)) 

target_test_surv <- DoSurvialAnalysis(target_test_score,
                  time = `Event Free Survival Time in Days`,
                  status = status,
                  predictor = score,
                  group_by = score_bin,
                  description = "TARGET Test Data",
                  lasso = target_model)

target_test_surv_os <- DoSurvialAnalysis(target_test_score,
                                      time = `Overall Survival Time in Days`,
                                      status = status,
                                      predictor = score,
                                      group_by = score_bin,
                                      description = "TARGET Test Data",
                                      lasso = target_model)

### Test with TCGA on TARGET model ----
tcga_data <- read_tsv(here("cibersort_in/tcga_cibersort.txt"))
tcga_data <- transposeDF(tcga_data, rowname_col = "case_submitter_id")

annotated_tcga <- read_tsv(here("data/tcga/clinical.cart.2021-04-22/clinical.tsv"), na = "'--") %>%
  right_join(tcga_data) %>%
  mutate(score = UseLASSOModelCoefs(cur_data(), coef(target_model)),
         score_bin = if_else(score >= median(score), "High", "Low"),
         status = if_else(vital_status == "Dead", 1, 0))

tcga_surv <- DoSurvialAnalysis(annotated_tcga,
                  time = days_to_death,
                  status = status,
                  predictor = score,
                  group_by = score_bin,
                  description = "TCGA Data (4.8% features present)",
                  lasso = target_model)

### Test with BeatAML on TARGET model ----

beat_aml <- read_tsv(here("data/aml_ohsu_2018/data_RNA_Seq_expression_cpm.txt"))

beat_aml_data <- beat_aml %>%
  dplyr::select(-Entrez_Gene_Id) %>%
  transposeDF(rowname_col = SAMPLE_ID) %>%
  mutate(across(.cols = !starts_with("SAMPLE_ID"), .fns = as.numeric)) %>%
  inner_join(beat_aml_clinical2) %>%
  filter(FLT3_ITD_CONSENSUS_CALL == "Positive") %>%
  mutate(across(.cols = !any_of(c("FLT3_ITD_CONSENSUS_CALL", "SAMPLE_ID")), .fns = as.numeric),
         OS_DAYS = OS_MONTHS * 365.25 / 12,
         score = UseLASSOModelCoefs(cur_data(), coef(target_model)),
         score_bin = if_else(score >= mean(score), "High", "Low"))

beataml_surv <- DoSurvialAnalysis(beat_aml_data,
                  time = OS_DAYS,
                  status = status,
                  predictor = score,
                  group_by = score_bin,
                  description = "Beat-AML Data (85.7% features present)",
                  lasso = target_model)

### Print Plots ----
pdf(here("plots/features_EFS_survival.pdf"))
target_train_surv$plot
target_train_surv_os$plot
target_test_surv$plot
target_test_surv_os$plot
tcga_surv$plot
beataml_surv$plot
graphics.off()
