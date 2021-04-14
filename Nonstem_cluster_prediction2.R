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
                       assays = c("ADT", "SCT"))

meta <- seurat[["patient_id"]] %>% rownames_to_column()
val <- read_excel(here("clinical_info/TARGET_AML_ClinicalData_Validation_20181213.xlsx"))
dis <- read_excel(here("clinical_info/TARGET_AML_ClinicalData_Discovery_20181213.xlsx"))
seq <- read_excel(file.path("../preprocessing/sample_info/Global Demultiplexing and Annotation.xlsx"))

UseLASSOModelCoefs <- function(x, coef) {
  mod2 <- coef[-1, ] %>% as.data.frame()
  w_vec <- list(0)
  mod2 <- mod2[which(mod2[1] != 0), ,drop = FALSE]
  
  for (i in seq_along(rownames(mod2))) {
    col <- rownames(mod2)[i]
    vec <- x[[col]]
    w <- mod2[i, 1]
    w_vec[[i]] <- vec * w
  }
  res <- w_vec %>% purrr::reduce(rbind) %>% colSums() + coef[1, ]
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
diagnosis <- readRDS(here("seurat_diagnosis_nonstem.Rds"))

FindClusterFreq(diagnosis, c("patient_id", "prognosis"), "clusters") %>%
  group_by(cluster) %>%
  summarise(pval = wilcox.test(freq ~ prognosis)$p.value, .groups = "keep") %>%
  arrange(pval)

clusters <- c("17", "3", "29", "5", "23", "18", "33", "26")

freq <- FindClusterFreq(diagnosis, c("patient_id", "prognosis"), "sct_clusters") %>%
  filter(cluster %in% clusters) %>%
  group_by(cluster) %>%
  dplyr::select(patient_id, freq, cluster, prognosis)

library(viridis)

pdf(file = here("plots/UMAP_and_freq_boxplot%01d.pdf"), onefile = FALSE, width = 10)
cells <- map(clusters, ~ WhichCells(diagnosis, idents = .x))
names(cells) <- clusters
UMAPPlot(diagnosis, split.by = "prognosis", cells.highlight = cells, cols.highlight = viridis(8)) 

ggplot(data = freq, mapping = aes(x = cluster, y = freq, fill = prognosis)) +
  geom_boxplot() +
  theme_classic() +
  xlab("Cluster") +
  ylab("Frequency per patient")
graphics.off()

summarise(freq,
          xbar = mean(freq),
          med = median(freq),
          q1 = quantile(freq)[2],
          q3 = quantile(freq)[4])

ChiSquarePValue <- function(survdiff) {
  df <- (sum(1 * (survdiff$exp > 0))) - 1
  pchisq(survdiff$chisq, df, lower.tail = FALSE)
}

survival_models <- freq %>%
  right_join(patient_info) %>%
  mutate(status = if_else(`First Event` == "Relapse", 1, 0),
                cluster_risk = if_else(freq >= mean(freq), "High", "Low")) %>%
  dplyr::select(`Event Free Survival Time in Days`, status, cluster_risk, patient_id, freq) %>%
  remove_missing() %>%
  nest() %>%
  mutate(survival = map(data, ~ survdiff(Surv(`Event Free Survival Time in Days`, status) ~ cluster_risk, data = .x))) %>%
  mutate(chi_sq = map_dbl(survival, ~ .x[["chisq"]]),
         freq = map_dbl(data, ~ mean(.x[["freq"]])),
         maxfreq = map_dbl(data, ~ max(.x[["freq"]])),
         minfreq = map_dbl(data, ~ min(.x[["freq"]])),
         p_val = map_dbl(survival, ChiSquarePValue)) %>%
  arrange(desc(chi_sq)) 


markers <- FindMarkers(diagnosis, ident.1 = survival_models$cluster[1:3])
markers_filtered <- markers[abs(markers$avg_log2FC) >= 0.5, ]

expression <- AverageExpression(diagnosis, group.by = "patient_id", features = rownames(markers))$SCT

clinical2 <- dplyr::select(clinical, patient_id, `Overall Survival Time in Days`,
                           `WBC at Diagnosis`, `Bone marrow leukemic blast percentage (%)`,
                           `Peripheral blasts (%)`, `Cytogenetic Complexity`,
                           `FLT3/ITD allelic ratio`, `Event Free Survival Time in Days`)

lasso_data <- expression %>%
  as.data.frame() %>%
  t() %>%
  as_tibble() %>%
  mutate(patient_id = colnames(expression), .before = 1) %>%
  inner_join(clinical2) 

lasso_data2 <- lasso_data[colnames(lasso_data) %in% colnames(target_data)]

vars <- dplyr::select(lasso_data2, -where(is.character), -`Event Free Survival Time in Days`, -`Overall Survival Time in Days`) %>% as.matrix()
pred <- dplyr::pull(lasso_data2, `Event Free Survival Time in Days`)
sc_model <- cv.glmnet(vars, pred)
sc_model_final <- glmnet(vars, pred, lambda = sc_model$lambda.min)

target <- read_tsv(file = here("target_data.txt"))

target_data <- data.table::transpose(target) %>%
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
  dplyr::select(-`FLT3/ITD positive?`, -patient_id) %>%
  mutate(across(.cols = !`First Event`, .fns = as.numeric)) 

sc_model_target_filtered <- coef(sc_model_final)[rownames(coef(sc_model_final)) %in% colnames(target_data), , drop = FALSE]
sc_model_target_filtered <- coef(sc_model_final)[1,,drop = FALSE] %>% rbind(sc_model_target_filtered)

target_test_score <- target_data %>%
  mutate(score = UseLASSOModelCoefs(.data, sc_model_target_filtered),
         score_bin = if_else(score >= mean(score, na.rm = TRUE), "High", "Low"),
         status = if_else(`First Event` == "relapse", 1, 0)) %>%
  dplyr::select(score, score_bin, `Event Free Survival Time in Days`, status)

score_survival <- survfit(Surv(`Event Free Survival Time in Days`, status) ~ score_bin, data = target_test_score)
survdiff(Surv(`Event Free Survival Time in Days`, status) ~ score_bin, data = target_test_score)

ggsurvplot(data = target_test_score,
           fit = score_survival, 
           xlab = "Days", 
           ylab = "Event Free Survival Probibility",
           risk.table = TRUE)

expression_t[colnames(expression_t) %in% colnames(target_data)]

pats <- colnames(target) %>%
  as_tibble() %>%
  separate(col = value, into = c(NA, NA, "patient_id", NA, NA)) %>%
  remove_missing() %>%
  pull(patient_id)

no_pats <- clinical[clinical$patient_id %!in% pats, ] %>% filter(`FLT3/ITD positive?` %in% c("YES", "Yes", "yes")) %>% pull(patient_id) %>% unique()
clinical[clinical$`FLT3/ITD positive?` %in% c("YES", "Yes", "yes"), ] %>% nrow()

sra <- read_tsv(here("SraRunTable.txt"))

sra %>%
  filter(Assay_Type == "RNA-Seq") %>%
  separate(Sample_Name, into = c(NA, NA, "patient_id", NA, NA)) %>%
  left_join(clinical) %>%
  filter(`FLT3/ITD positive?` %in% c("YES", "Yes", "yes") & clinical$patient_id %!in% pats) %>%
  pull()
