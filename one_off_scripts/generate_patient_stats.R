library(tidyverse)
library(here)
library(broom)

source(here("cli/lib/prepare_clin_info.R"))

clinical # TARGET
tcga_ann2 # TCGA
beat_aml_clinical2 # BeatAML

no <- c("no", "NO", "No", "Unknown", "unknown")

target_summary <- clinical %>%
  mutate(across(where(is.character), str_to_lower)) %>%
  mutate(
    across(where(is.character), as.factor),
    subset = if_else(
      .data[["FLT3/ITD positive?"]] %in% no,
      if_else(
        .data[["CEBPA mutation"]] %in% no,
        "NEG",
        "CEBPA"
      ),
      "FLT3"
    ),
    `FLT3/ITD allelic ratio` = as.numeric(`FLT3/ITD allelic ratio`)
  ) %>%
  group_by(subset) %>%
  select(-Comment) %>%
  distinct() %>%
  nest() %>%
  mutate(
    summary = map(data, summary)
  )

target_summary$summary

tcga_ann2 %>%
  filter(flt3_status == "Positive") %>%
  mutate(across(where(is.character), as.factor)) %>%
  summary()

beat_aml_clinical %>%
  filter(FLT3_ITD_CONSENSUS_CALL == "Positive") %>%
  mutate(across(where(is.character), as.factor)) %>%
  summary()

library(Seurat)
library(SeuratDisk)

cebpa_patients <- target_summary$data[[3]] %>%
  pull(USI) %>%
  str_to_upper()

cebpa_patients <- c("PAUUWT", "PAUVMN", "PAVCNG", "PAWAMM", "PAWHML", "PAWUWJ", "PAWZZN")



seurat <- LoadH5Seurat("data/preprocessed.h5Seurat", assays = FALSE)

seurat[[]] %>%
  filter(timepoint == "Diagnosis") %>%
  group_by(prognosis) %>%
  summarise(n = length(unique(patient_id)))

seurat[[]] %>%
  filter(timepoint == "Diagnosis") %>%
  mutate(CEBPA = if_else(patient_id %in% cebpa_patients, "CEBPA", "FLT3")) %>%
  group_by(patient_id, stemness, CEBPA) %>%
  summarise(
    cells = n(),
    percent_mt = median(percent.mt),
    across(starts_with("n"), median)
  ) %>%
  select(-contains(c("AML", "PAN", "ADT"))) %>%
  write_tsv("one_off_scripts/sequencing_info.tsv")

seurat[[]] %>%
  filter(timepoint == "Diagnosis") %>%
  mutate(CEBPA = if_else(patient_id %in% cebpa_patients, "CEBPA", "FLT3")) %>%
  group_by(prognosis) %>%
  summarise(
    n = n()
  )

our_pats <- seurat[[]] %>%
  filter(timepoint == "Diagnosis") %>%
  pull(patient_id) %>%
  unique()

read_excel("clinical_info/AAML19B3Q_data_transfer.xlsx", sheet = 2) %>%
  filter(USI %in% our_pats) %>%
  mutate(
    across(where(is.character), as.factor),
    `Allelic ratio` = as.numeric(`Allelic ratio`),
    `Time to relapse from study entry in days` = as.numeric(as.character(`Time to relapse from study entry in days`))
  ) %>%
  distinct() %>%
  summary()
