library(survminer)
library(survival)
library(here)
library(broom)
library(tidyverse)

source(here("cli/lib/prepare_clin_info.R"))

flt3 <- clinical %>%
  filter(`FLT3/ITD positive?` == "Yes") %>%
  mutate(
    event = if_else(`First Event` == "Censored", 0, 1),
    AR    = if_else(`FLT3/ITD allelic ratio` > 0.4, "AR+", "AR-")
  ) %>%
  sjlabelled::set_na(na = "Unknown")

surv_form <- list(
  AR    = Surv(`Event Free Survival Time in Days`, event) ~ AR,
  PM    = Surv(`Event Free Survival Time in Days`, event) ~ `FLT3 PM`,
  MLL   = Surv(`Event Free Survival Time in Days`, event) ~ MLL,
  SCT   = Surv(`Event Free Survival Time in Days`, event) ~ `SCT in 1st CR`,
  NPM   = Surv(`Event Free Survival Time in Days`, event) ~ `NPM mutation`,
  CEBPA = Surv(`Event Free Survival Time in Days`, event) ~ `CEBPA mutation`,
  CK8   = Surv(`Event Free Survival Time in Days`, event) ~ `c-Kit Mutation Exon 8`,
  CK17  = Surv(`Event Free Survival Time in Days`, event) ~ `c-Kit Mutation Exon 17`
)

surv_form <- Surv(`Event Free Survival Time in Days`, event) ~
AR + MLL + `NPM mutation` + `WT1 mutation` +
  `CEBPA mutation` + `c-Kit Mutation Exon 8` + `c-Kit Mutation Exon 17`

fit <- coxph(surv_form, data = flt3)
plot1 <- ggforest(fit)

surv_form <- Surv(`Event Free Survival Time in Days`, event) ~
AR + MLL + `WT1 mutation` +
  `CEBPA mutation` + `c-Kit Mutation Exon 8` + `c-Kit Mutation Exon 17`

flt3 %>%
  filter(!(`NPM mutation` == "Yes" & AR == "AR-")) %>%
  coxph(surv_form, .)

surv_form <- Surv(`Event Free Survival Time in Days`, event) ~ itd

all <- clinical %>%
  mutate(
    event = if_else(`First Event` == "Censored", 0, 1),
    itd = str_to_lower(`FLT3/ITD positive?`)
  ) %>%
  sjlabelled::set_na(na = "unknown")

fit <- coxph(surv_form, all)

fit2efs <- survfit(Surv(`Event Free Survival Time in Days`, event) ~ itd, all)
fit2os <- survfit(Surv(`Overall Survival Time in Days`, event) ~ itd, all)
plot2 <- ggsurvplot(
  list(fit2efs, fit2os),
  data = all,
  conf.int = TRUE,
  palette = c("#E7B800", "#2E9FDF"),
  risk.table = "abs_pct",
  surv.median.line = "hv"
)

pdf("one_off_scripts/std_bio_ggforest.pdf")
plot1
plot2
graphics.off()
