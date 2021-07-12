library(here)
library(readxl)
library(xlsx)
library(tidyverse)

source(here("lib/functions.R"))
### TARGET ----

val <- read_excel(
  here("clinical_info/TARGET_AML_ClinicalData_Validation_20181213.xlsx"),
  col_types = cols()
)
dis <- read_excel(
  here("clinical_info/TARGET_AML_ClinicalData_Discovery_20181213.xlsx"),
  col_types = cols()
)
cog <- read.xlsx(here("clinical_info/AAML19B3Q_data_transfer.xlsx"),
  sheetIndex = 2,
  password = "AAML19B3Q"
) %>%
  sjlabelled::set_na(na = ".") %>%
  as_tibble() %>%
  rename(
    init_treatment_arm = Treatment.arm.at.enrollment,
    fin_treatment_arm = AAML1031..Final.treatment.arm.assignment,
    `Overall Survival Time in Days` = Days.to.OS.from.study.entry,
    `FLT3/ITD allelic ratio` = Allelic.ratio,
    `WBC at Diagnosis` = WBC..x10.3.MicroLiter..,
    `Event Free Survival Time in Days` = Time.to.relapse.from.study.entry.in.days
  ) %>%
  mutate(
    `FLT3/ITD positive?` = if_else(
      FLT3.results == "Internal tandem duplication",
      "Yes",
      "No"
    ),
    `CEBPA mutation` = if_else(CEBPA.mutation.status == "Positive",
      "Yes",
      "No"
    ),
    `NPM mutation` = if_else(Necleophosmin..NPM..mutation.status == "Positive",
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
    "../preprocessing/sample_info/Global Demultiplexing and Annotation.xlsx"
  )
)

clinical <- bind_rows(val, dis) %>%
  distinct() %>%
  separate(`TARGET USI`, into = c(NA, NA, "USI"), sep = "-")

### TCGA ----
tcga_ann <- read_tsv(
  here("data/tcga/clinical.cart.2021-04-22/clinical.tsv"),
  na = "'--"
)

tcga_ann2 <- read_tsv(
  here("clinical_info/nationwidechildrens.org_clinical_patient_laml.txt"),
  skip = 1,
  na = c("[Not Available]", "[Not Applicable]")
) %>%
  filter(bcr_patient_uuid != "CDE_ID:") %>%
  map_df(parse_guess) %>%
  mutate(
    flt3_status = str_extract(
      molecular_analysis_abnormality_testing_result,
      "FLT3 Mutation [:alpha:]+"
    )
  ) %>%
  separate(flt3_status, into = c(NA, NA, "flt3_status"), sep = " ") %>%
  rename(`Case ID` = bcr_patient_barcode) %>%
  left_join(tcga_ann) %>%
  mutate(case_submitter_id = `Case ID`)

### BeatAML ----
beat_aml_clinical <- read_tsv(
  here("data/aml_ohsu_2018/data_clinical_sample.txt"),
  skip = 4
)
beat_aml_survival <- read_tsv(
  here("data/aml_ohsu_2018/KM_Plot__Overall_Survival__(months).txt")
)

beat_aml_clinical2 <- beat_aml_clinical %>%
  inner_join(beat_aml_survival, by = c("PATIENT_ID" = "Patient ID")) %>%
  mutate(
    status = if_else(SAMPLE_TIMEPOINT == "Relapse", 1, 0),
    OS_DAYS = OS_MONTHS * (365.25 / 12)
  ) %>%
  select(
    SAMPLE_ID,
    FLT3_ITD_CONSENSUS_CALL,
    OS_DAYS,
    PB_BLAST_PERCENTAGE,
    status
  ) %>%
  rename(`Peripheral blasts (%)` = PB_BLAST_PERCENTAGE)
