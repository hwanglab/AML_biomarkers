library(here)
library(readxl)
library(tidyverse)

val <- read_excel(
  here("clinical_info/TARGET_AML_ClinicalData_Validation_20181213.xlsx")
)
dis <- read_excel(
  here("clinical_info/TARGET_AML_ClinicalData_Discovery_20181213.xlsx")
)
cog <- read_xlsx(here("clinical_info/AAML19B3Q_data_transfer.xlsx"),
  sheet = 2,
) %>%
  sjlabelled::set_na(na = ".") %>%
  as_tibble() %>%
  rename(
    init_treatment_arm = `Treatment arm at enrollment`,
    fin_treatment_arm = `AAML1031: Final treatment arm assignment`,
    `Overall Survival Time in Days` = `Days to OS from study entry`,
    `FLT3/ITD allelic ratio` = `Allelic ratio`,
    `WBC at Diagnosis` = `WBC (x10^3/MicroLiter)`,
    `Event Free Survival Time in Days` = `Time to relapse from study entry in days`
  ) %>%
  mutate(
    `FLT3/ITD positive?` = if_else(
      `FLT3 results` == "Internal tandem duplication",
      "Yes",
      "No"
    ),
    `CEBPA mutation` = if_else(`CEBPA mutation status` == "Positive",
      "Yes",
      "No"
    ),
    `NPM mutation` = if_else(`Necleophosmin (NPM) mutation status` == "Positive",
      "Yes",
      "No"
    ),
    `Event Free Survival Time in Days` = as.numeric(
      `Event Free Survival Time in Days`
    ),
    `Event Free Survival Time in Days` = if_else(
      `Event Free Survival Time in Days` == ".",
      5 * 365.25,
      `Event Free Survival Time in Days`
    )
  ) %>%
  select(
    init_treatment_arm, fin_treatment_arm, `Overall Survival Time in Days`,
    `FLT3/ITD positive?`, `WBC at Diagnosis`, `FLT3/ITD positive?`,
    `NPM mutation`, `CEBPA mutation`, USI
  )

seq <- read_excel(
  file.path(
    stringr::str_replace(Sys.getenv("AML_DATA"), "outs", "sample_info"),
    "Global Demultiplexing and Annotation.xlsx"
  )
)

clinical <- bind_rows(val, dis) %>%
  distinct() %>%
  separate(`TARGET USI`, into = c(NA, NA, "USI"), sep = "-")
