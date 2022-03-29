library(survminer)
library(survival)
library(here)
library(broom)
library(tidyverse)

source(here("cli/lib/prepare_clin_info.R"))

flt3 <- clinical %>% filter(`FLT3/ITD positive?` == "Yes") %>%
    mutate(
        event = if_else(`First Event` == "Relapse",       1,      0),
        AR    = if_else(`FLT3/ITD allelic ratio`   > 0.2, "AR+",  "AR-"),
        MRD   = if_else(`MRD % at end of course 1` > 0.4, "MRD+", "MRD-")
    )

surv_form <- list(
    AR  = Surv(`Event Free Survival Time in Days`, event) ~ AR,
    PM  = Surv(`Event Free Survival Time in Days`, event) ~ `FLT3 PM`,
    MRD = Surv(`Event Free Survival Time in Days`, event) ~ MRD
)

map(surv_form, ~coxph(.x, data = flt3)) %>% map(summary)

surv_form <- Surv(`Event Free Survival Time in Days`, event) ~ MRD + `FLT3 PM` + AR
coxph(surv_form, data = flt3) %>% summary()
fit <- coxph(surv_form, data = flt3)
