## Figure 2 of Paper

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

#signature <- read_table("data/lsc_de_ng_et_al.txt", col_names = FALSE)[[1]]

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
  as_tibble() %>%
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

pdf("plots/LSC48_prediction.pdf")
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

## Test on TCGA ----

tcga_data <- read_tsv(here("cibersort_in/tcga_data.txt"))
tcga_data <- transposeDF(tcga_data, rowname_col = "case_submitter_id")

annotated_tcga <- read_tsv(
  here("data/tcga/clinical.cart.2021-04-22/clinical.tsv"),
  na = "'--"
) %>%
  right_join(tcga_data) %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
    score_bin = if_else(score >= median(score), "High", "Low"),
    status = if_else(vital_status == "Dead", 1, 0)
  )

tcga_surv <- DoSurvivalAnalysis(
  annotated_tcga,
  days_to_death,
  status,
  score,
  group_by = score_bin,
  description = "TCGA (33% features present)",
  lasso = lasso_model
)

## Test on BeatAML ----
beat_aml <- read_tsv(here("data/aml_ohsu_2018/data_RNA_Seq_expression_cpm.txt"))

beat_aml_clinical <- read_tsv(
  here("data/aml_ohsu_2018/data_clinical_sample.txt"),
  skip = 4
)
beat_aml_survival <- read_tsv(
  here("data/aml_ohsu_2018/KM_Plot__Overall_Survival__(months).txt")
)
beat_aml_clinical2 <- inner_join(
  beat_aml_clinical,
  beat_aml_survival,
  by = c("PATIENT_ID" = "Patient ID")
) %>%
  mutate(status = if_else(SAMPLE_TIMEPOINT == "Relapse", 1, 0)) %>%
  dplyr::select(
    SAMPLE_ID,
    FLT3_ITD_CONSENSUS_CALL,
    OS_MONTHS,
    PB_BLAST_PERCENTAGE,
    status
  ) %>%
  rename(`Peripheral blasts (%)` = PB_BLAST_PERCENTAGE)


beat_aml_data <- beat_aml %>%
  dplyr::select(-Entrez_Gene_Id) %>%
  transposeDF(rowname_col = SAMPLE_ID) %>%
  mutate(across(.cols = !starts_with("SAMPLE_ID"), .fns = as.numeric)) %>%
  inner_join(beat_aml_clinical2) %>%
  filter(FLT3_ITD_CONSENSUS_CALL == "Positive") %>%
  mutate(
    across(
      .cols = !any_of(c("FLT3_ITD_CONSENSUS_CALL", "SAMPLE_ID")),
      .fns = as.numeric
    ),
    OS_DAYS = OS_MONTHS * 365.25 / 12,
    score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
    score_bin = if_else(score >= mean(score), "High", "Low")
  )

beat_surv <- DoSurvivalAnalysis(
  beat_aml_data,
  OS_DAYS,
  status,
  score,
  group_by = score_bin,
  description = "BeatAML (67% features present)",
  lasso = lasso_model
)

# print plots
pdf(file = here("plots/CIBERSORTless_TARGET_trained.pdf"))
target_train$plot
target_test$plot
target_train_os$plot
target_test_os$plot
tcga_surv$plot
beat_surv$plot
graphics.off()

stats <- list(
  "TARGET (Train)" = target_train,
  "TARGET (Test)" = target_test,
  "TARGET (Train) OS" = target_train_os,
  "TARGET (Test) OS" = target_test_os,
  "TCGA" = tcga_surv,
  "BeatAML" = beat_surv
)

stat_res <- map(stats, ~ .x[["stat"]]) %>%
  map_dfr(broom::tidy) %>%
  mutate(dataset = names(stats))

write_tsv(
  stat_res,
  here("outs/TARGET_LASSO_model_survival_without_CIBERSORTx.tsv")
)