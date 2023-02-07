
library(Seurat)
library(tidyverse)

# Load the data
seurat <- readRDS("outs/flt3_cebpa/cache/seurat_d307f778811553cd7be8a49668d726e4.rds")

markers <- FindAllMarkers(seurat, only.pos = TRUE, test.use = "LR")

# Write the results to a file
write_tsv(markers, "one_off_scripts/cluster_marker_model_markers_nonstem.tsv")

markers <- read_tsv("one_off_scripts/cluster_marker_model_markers_nonstem.tsv")

# use chi_square clusters
signature <- markers %>%
    filter(cluster %in% c(1, 4, 5, 8, 13, 16, 17, 18, 19, 21, 22, 24, 25, 27, 29, 30, 33, 35, 36)) %>%
    group_by(cluster) %>%
    slice(1:10) %>% 
    pull(gene)

library(here)
library(readxl)
library(survival)
library(survminer)
library(scorecard)
library(glmnet)
library(data.table)
library(tidyverse)

source(here("cli/lib/lasso.R"))
source(here("cli/lib/utils.R"))

target <- read_tsv(file = here("cibersort_in/target_data.txt"))
val <- read_excel(here("clinical_info/TARGET_AML_ClinicalData_Validation_20181213.xlsx"))
dis <- read_excel(here("clinical_info/TARGET_AML_ClinicalData_Discovery_20181213.xlsx"))


clinical <- bind_rows(val, dis) %>%
  distinct() %>%
  separate(`TARGET USI`, into = c(NA, NA, "patient_id"), sep = "-")

clinical2 <- dplyr::select(
  clinical, `WBC at Diagnosis`, `Bone marrow leukemic blast percentage (%)`,
  `Peripheral blasts (%)`, `Cytogenetic Complexity`,
  `FLT3/ITD allelic ratio`, `Event Free Survival Time in Days`,
  `FLT3/ITD positive?`, patient_id, `First Event`,
  `Overall Survival Time in Days`
)

rn <- target[[1]]
cn <- colnames(target)

target <- target[, -1]

rownames(target) <- colnames(target) <- NULL

target_before_split <- target %>%
  as.data.table() %>%
  transpose() %>%
  as_tibble(.name_repair = "minimal") %>%
  set_names(rn) %>%
  mutate(
    across(.fns = as.numeric),
    samples = cn[-1]
  ) %>%
  separate(col = "samples", into = c(NA, NA, "patient_id", NA, NA)) %>%
  inner_join(clinical2) %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "yes") %>%
  dplyr::select(-`FLT3/ITD positive?`, -patient_id) %>%
  mutate(across(.cols = !`First Event`, .fns = as.numeric)) %>%
  remove_missing()

target_split <- split_df(
  target_before_split,
  y = "Event Free Survival Time in Days",
  ratios = c(0.7, 0.3)
)

## Train on TARGET (train) ----
training <- target_split$train %>%
  select(any_of(signature), `Event Free Survival Time in Days`)
lasso_model <- DoGLMAlpha(
  training, "Event Free Survival Time in Days", 1, 10,
  k = 1000, save = FALSE
)[[1]]

## Test on TARGET (test) ----
test <- target_split$test %>%
  as_tibble() %>%
  mutate(
    status = if_else(`First Event` != "censored", 1, 0),
    `First Event` = NULL
  ) %>%
  mutate(
    score = predict(lasso_model, cur_data()),
    score_bin = if_else(score >= median(score), "High", "Low")
  )

fit <- survfit(Surv(`Event Free Survival Time in Days`, status) ~ score_bin, data = test)

coxph(Surv(`Event Free Survival Time in Days`, status) ~ score_bin, data = test)

pdf("plots/top_n_genes_prediction.pdf")
survminer::ggsurvplot(
  fit,
  ylab = "Survival Probability",
  pval = T,
  xlab = "Survival Time",
  font.subtitle = 8,
  risk.table = "abs_pct",
  palette = c("#E06C9F", "#23395B"),
  conf.int = TRUE
)
graphics.off()
