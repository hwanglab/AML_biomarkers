library(survminer)
library(survival)
library(here)
library(broom)
library(tidyverse)

source(here("cli/lib/prepare_clin_info.R"))

flt3 <- clinical %>%
  mutate(
    event = if_else(`First Event` == "Censored", 0, 1),
    flt3_pos = str_to_lower(`FLT3/ITD positive?`),
  ) %>%
  sjlabelled::set_na(na = "Unknown") %>%
  filter(flt3_pos != "unknown") %>%
  rename(
    efs = `Event Free Survival Time in Days`,
    age = `Age at Diagnosis in Days`,
    os = `Overall Survival Time in Days`
  ) %>%
  mutate(
    efs = efs / 365.25 * 12,
    os = os / 365.25 * 12
  )

surv_form <- Surv(efs, event) ~ flt3_pos + Gender + Race + Ethnicity + age

fit <- coxph(surv_form, data = flt3)
plot1 <- ggforest(fit)

plot2 <- ggadjustedcurves(fit, variable = "flt3_pos", data = as.data.frame(flt3))

surv_form <- Surv(os, event) ~ flt3_pos + Gender + Race + Ethnicity + age

fit <- coxph(surv_form, data = flt3)
plot3 <- ggforest(fit)

plot4 <- ggadjustedcurves(fit, variable = "flt3_pos", data = as.data.frame(flt3))
pdf("one_off_scripts/flt3_itd.pdf")
plot1
plot2
plot3
plot4
graphics.off()
