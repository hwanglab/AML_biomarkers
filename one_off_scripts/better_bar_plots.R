library(Seurat)
library(tidyverse)
library(viridis)
library(ggsignif)
library(survival)

PrepareData <- function(seurat) {
  seurat[[c("seurat_clusters", "patient_id", "prognosis")]] %>%
    group_by(seurat_clusters, patient_id, prognosis) %>%
    summarise(n = n()) %>%
    group_by(patient_id) %>%
    mutate(f = n / sum(n) * 100)
}

PrepareData2 <- function(seurat) {
  seurat[[c("seurat_clusters", "patient_id", "prognosis")]] %>%
    group_by(seurat_clusters, patient_id, prognosis) %>%
    summarise(n = n()) %>%
    group_by(seurat_clusters) %>%
    mutate(f = n / sum(n) * 100)
}


PrepareSummaryStatistics <- function(df) {
  df %>%
  group_by(seurat_clusters, prognosis) %>%
  summarise(
    sd = sd(f, na.rm = TRUE),
    count = n(),
    mean = mean(f),
    median = median(f),
    sem = sd / sqrt(count)
  )
}

PrintPlot <- function(df, stemness, data) {
  ggplot(data = df, mapping = aes(x = prognosis, fill = prognosis, y = mean)) +
  geom_col() +
  facet_wrap(~seurat_clusters, scales = "free") +
  theme_classic() +
  geom_errorbar(
    aes(ymin = mean - sem, ymax = mean + sem),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  geom_jitter(data = data, mapping = aes(x = prognosis, y = f)) +
  labs(x = NULL, caption = stemness) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

seurat <- list(
  stem = readRDS("outs/flt3_cebpa_stem/cache/seurat_365ef0e0d8dbade1c5eb33690b90ab6a.rds"),
  nonstem = readRDS("outs/flt3_cebpa/cache/seurat_d307f778811553cd7be8a49668d726e4.rds")
)

data <- map(seurat, PrepareData)
data2 <- map(seurat, PrepareData2)

stats <- map(data, PrepareSummaryStatistics)

plots <- pmap(
  tibble(stats = stats, stemness = names(stats), data = data),
  ~PrintPlot(..1, ..2, ..3)
)

idx <- "stem"
pdf("one_off_scripts/bar_plot1.pdf")
ggplot(data = stats[[idx]], mapping = aes(x = prognosis, fill = prognosis, y = mean)) +
  geom_col() +
  facet_wrap(~seurat_clusters, scales = "free") +
  theme_classic() +
  geom_errorbar(
    aes(ymin = mean - sem, ymax = mean + sem),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  geom_jitter(data = data[[idx]], mapping = aes(x = prognosis, y = f)) +
  labs(x = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

data[[idx]] %>% 
  group_by(seurat_clusters) %>% 
  mutate(f = n / sum(n) * 100) %>%
  ungroup() %>%
  arrange(prognosis) %>%
  mutate(patient_id = factor(patient_id, levels = unique(patient_id))) %>%
  ggplot(mapping = aes(x = seurat_clusters, y = f, fill = patient_id)) +
  geom_col() +
  scale_alpha_manual(values = c("Poor" = 0.5, "Good" = 1)) +
  theme_classic()
graphics.off()

idx <- "nonstem"
pdf("one_off_scripts/bar_plot2.pdf")
ggplot(data = stats[[idx]], mapping = aes(x = prognosis, fill = prognosis, y = mean)) +
  geom_col() +
  facet_wrap(~seurat_clusters, scales = "free") +
  theme_classic() +
  geom_errorbar(
    aes(ymin = mean - sem, ymax = mean + sem),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  geom_jitter(data = data[[idx]], mapping = aes(x = prognosis, y = f)) +
  labs(x = NULL) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

data[[idx]] %>% 
  group_by(seurat_clusters) %>% 
  mutate(f = n / sum(n) * 100) %>%
  ungroup() %>%
  arrange(prognosis) %>%
  mutate(patient_id = factor(patient_id, levels = unique(patient_id))) %>%
  ggplot(mapping = aes(x = seurat_clusters, y = f, fill = patient_id)) +
  geom_col() +
  scale_alpha_manual(values = c("Poor" = 0.5, "Good" = 1)) +
  theme_classic()
graphics.off()

t_test <- purrr::possibly(wilcox.test, otherwise = NA_real_)
survdiff_safe <- purrr::possibly(survdiff, otherwise = NA_real_)
pchisq_safe <- purrr::possibly(pchisq, otherwise = NA_real_)
cox_safe <- purrr::possibly(coxph, otherwise = NA_real_)
tidy_safe <- purrr::possibly(broom::tidy, otherwise = NA_real_)

source("cli/lib/prepare_clin_info.R")

clinical <- select(cog, USI, efs = cog_efs)

tibble(set = names(data), info = data) %>%
  group_by(set) %>%
  unnest(info) %>%
  group_by(set, seurat_clusters) %>%
  left_join(clinical, by = c("patient_id" = "USI")) %>%
  filter(!is.na(efs)) %>%
  mutate(event = if_else(prognosis == "Poor", 1, 0)) %>%
  nest() %>%
  mutate(
    # Wilcoxon rank-sum test
    t              = map(data, ~t_test(f ~ prognosis, data = .x)),
    t_test_p_value = map_dbl(t, pluck, "p.value", .default = NA_real_),

    # Linear regression
    lm        = map(data, ~ lm(f ~ efs, data = .x)),
    summary   = map(lm, summary),
    tidy      = map(summary, broom::tidy),
    r_sqaured = map_dbl(summary, "r.squared"),

    # Survival analysis
    surv         = map(
                    data, 
                    ~ survdiff_safe(Surv(efs, event) ~ f, data = .x)
    ),
    surv_p_value = map_dbl(
                    surv, 
                    ~pchisq_safe(.x$chisq, length(.x$n) - 1, lower.tail = FALSE)
    ),

    # Cox regression
    cox = map(data, ~ cox_safe(Surv(efs, event) ~ f, data = .x)),
    cox_p = map(cox, ~ tidy_safe(.x)),
    cox_p = map_dbl(cox_p, ~ pluck(.x, "p.value", .default = NA_real_))
  ) %>%
  ungroup() %>%
  hoist(tidy, "p.value", "term") %>%
  rename(lm_p_value = p.value) %>%
  mutate(
    lm_p_value = map2(lm_p_value, term, set_names),
    lm_n = map_dbl(lm, ~ ncol(table(model.frame(.x)))),
  ) %>%
  unnest_wider(lm_p_value) %>%
  mutate(
    efs = if_else(lm_n <= 3, NA_real_, efs),
    r_sqaured = if_else(lm_n <= 3, NA_real_, r_sqaured),
    n = map_int(data, nrow),
    n_rel = map_dbl(data, ~ nrow(filter(.x, prognosis == "Poor"))),
    n_fav = map_dbl(data, ~ nrow(filter(.x, prognosis == "Favorable"))),
  ) %>%
  select(-data, -lm, -summary, -tidy, -term, -t, -`(Intercept)`, -surv, -efs, -r_sqaured, -cox) %>%
  write_tsv("one_off_scripts/summary_stats_bar_plot.tsv")

tibble(set = names(data), info = data) %>%
  unnest_wider(info) %>%
  unnest(c(seurat_clusters, patient_id, prognosis, n, f)) %>%
  group_by(set, seurat_clusters) %>% 
  pivot_wider(names_from = prognosis, values_from = n) %>%
  summarise(
    across(c(Poor, Favorable), sum, na.rm = TRUE)
  ) %>%
  write_tsv("one_off_scripts/cells_per_cluster.tsv")

tibble(set = names(data2), info = data2) %>%
  unnest_wider(info) %>%
  unnest(c(seurat_clusters, patient_id, prognosis, n, f)) %>%
  group_by(set, seurat_clusters) %>%
  select(-n, -prognosis) %>%
  pivot_wider(names_from = patient_id, values_from = f) %>%
  mutate(across(where(is_double), ~ if_else(is.na(.x), 0, .x))) %>%
  write_tsv("one_off_scripts/cluster_by_patient_freq.tsv")
