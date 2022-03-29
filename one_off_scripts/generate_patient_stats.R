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
            "CEBPA"),
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
