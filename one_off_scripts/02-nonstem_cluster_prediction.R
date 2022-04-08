library(Seurat)
library(SeuratDisk)
library(here)
library(rsinglecell)
library(readxl)
library(xlsx)
library(biomaRt)
library(granulator)
library(survival)
library(survminer)
library(scorecard)
library(glmnet)
library(xfun)
library(msigdbr)
library(furrr)
library(EnhancedVolcano)
library(plotly)
library(tidyverse)

source(here("lib/functions.R"))

## Globals ----
sig_level <- 0.10
rerun_cellphone <- FALSE

## Prepare Clinical Data ----
### TARGET ----
val <- read_excel(here("clinical_info/TARGET_AML_ClinicalData_Validation_20181213.xlsx"))
dis <- read_excel(here("clinical_info/TARGET_AML_ClinicalData_Discovery_20181213.xlsx"))
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
    `FLT3/ITD positive?` = if_else(FLT3.results == "Internal tandem duplication",
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
    `Event Free Survival Time in Days` = as.numeric(`Event Free Survival Time in Days`),
    `Event Free Survival Time in Days` = if_else(`Event Free Survival Time in Days` == ".",
      5 * 365.25,
      `Event Free Survival Time in Days`
    )
  ) %>%
  select(
    init_treatment_arm, fin_treatment_arm, `Overall Survival Time in Days`,
    `FLT3/ITD positive?`, `WBC at Diagnosis`, `FLT3/ITD positive?`,
    `NPM mutation`, `CEBPA mutation`, USI
  )

seq <- read_excel(file.path("../preprocessing/sample_info/Global Demultiplexing and Annotation.xlsx"))

clinical <- bind_rows(val, dis) %>%
  distinct() %>%
  separate(`TARGET USI`, into = c(NA, NA, "USI"), sep = "-")

### TCGA ----
tcga_ann <- read_tsv(here("data/tcga/clinical.cart.2021-04-22/clinical.tsv"), na = "'--")

tcga_ann2 <- read_tsv(here("clinical_info/nationwidechildrens.org_clinical_patient_laml.txt"), skip = 1, na = c("[Not Available]", "[Not Applicable]")) %>%
  filter(bcr_patient_uuid != "CDE_ID:") %>%
  map_df(parse_guess) %>%
  mutate(flt3_status = str_extract(molecular_analysis_abnormality_testing_result, "FLT3 Mutation [:alpha:]+")) %>%
  separate(flt3_status, into = c(NA, NA, "flt3_status"), sep = " ") %>%
  rename(`Case ID` = bcr_patient_barcode) %>%
  left_join(tcga_ann) %>%
  mutate(case_submitter_id = `Case ID`)

### BeatAML ----
beat_aml_clinical <- read_tsv(here("data/aml_ohsu_2018/data_clinical_sample.txt"), skip = 4)
beat_aml_survival <- read_tsv(here("data/aml_ohsu_2018/KM_Plot__Overall_Survival__(months).txt"))
beat_aml_clinical2 <- beat_aml_clinical %>%
  inner_join(beat_aml_survival, by = c("PATIENT_ID" = "Patient ID")) %>%
  mutate(status = if_else(SAMPLE_TIMEPOINT == "Relapse", 1, 0)) %>%
  select(SAMPLE_ID, FLT3_ITD_CONSENSUS_CALL, OS_MONTHS, PB_BLAST_PERCENTAGE, status) %>%
  rename(`Peripheral blasts (%)` = PB_BLAST_PERCENTAGE)

## Load Bulk Data ----
target_data <- read_tsv(here("cibersort_in/target_data.txt")) %>%
  column_to_rownames(var = "GeneSymbol") %>%
  as.matrix()

tcga_data <- read_tsv(here("cibersort_in/tcga_data.txt")) %>%
  column_to_rownames(var = "GeneSymbol") %>%
  as.matrix()

beatAML_data <- read_tsv(here("cibersort_in/beat_aml.txt")) %>%
  column_to_rownames(var = "GeneSymbol") %>%
  as.matrix()

## Subset + Dimension Reductions ----
diagnosis <- cache_rds(
  expr = {
    seurat <- LoadH5Seurat(file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat"),
      assays = c("SCT", "RNA")
    )

    DefaultAssay(seurat) <- "SCT"

    meta <- seurat[["patient_id"]] %>% rownames_to_column()
    data <- filter(clinical, patient_id %in% pull(seq, patient_id))

    seurat@misc[["patient_data"]] <- data
    seurat@misc[["target_data"]] <- clinical

    diagnosis <- subset(seurat, timepoint == "Diagnosis" & stemness == "Nonstem")
    diagnosis <- DoDimensionReductions(diagnosis, batch_vars = c("seq_batch", "sort_batch"))

    diagnosis[["clusters"]] <- Idents(diagnosis)
    diagnosis
  },
  file = "02-seurat_diagnosis_nonstem.rds",
  hash = list(file.info(file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat")))
)

## Find DE Clusters ----
wilcox_clusters <- FindClusterFreq(diagnosis[[]], c("patient_id", "prognosis"), "clusters") %>%
  group_by(cluster) %>%
  summarise(pval = wilcox.test(freq ~ prognosis)$p.value, .groups = "keep") %>%
  filter(pval <= sig_level) %>%
  pull(cluster)

freq <- FindClusterFreq(diagnosis[[]], c("patient_id", "prognosis"), "clusters") %>%
  group_by(cluster) %>%
  dplyr::select(patient_id, freq, cluster)

summarise(freq,
  xbar = mean(freq),
  med = median(freq),
  q1 = quantile(freq)[2],
  q3 = quantile(freq)[4]
)

survival_models <- freq %>%
  right_join(clinical, by = c("patient_id" = "USI")) %>%
  mutate(
    status = if_else(`First Event` == "Relapse", 1, 0),
    cluster_risk = if_else(freq >= mean(freq), "High", "Low")
  ) %>%
  dplyr::select(`Event Free Survival Time in Days`, status, cluster_risk, patient_id, freq) %>%
  remove_missing() %>%
  filter(!near(freq, mean(freq))) %>%
  nest() %>%
  mutate(survival = map(data, ~ survdiff(Surv(`Event Free Survival Time in Days`, status) ~ cluster_risk, data = .x))) %>%
  mutate(
    chi_sq = map_dbl(survival, ~ .x[["chisq"]]),
    p_val = map_dbl(survival, CalculatePValues.chi),
    log_p = -log(p_val),
    freq = map_dbl(data, ~ mean(.x[["freq"]])),
    wilcox = if_else(cluster %in% wilcox_clusters, "wilcox_sig", "not_wilcox_sig"),
    chi_sig = if_else(p_val <= sig_level, "chi_sig", "not_chi_sig"),
    sd = map_dbl(data, ~ sd(.x[["freq"]])),
    mean_diff = map_dbl(data, ReturnDifferences)
  ) %>%
  arrange(p_val)

pdf(file = here("plots/sc_cluster_survival_analysis.pdf"), width = 3, height = 5)
ggplot(data = survival_models, mapping = aes(x = mean_diff, y = log_p)) +
  geom_hline(yintercept = -log(sig_level)) +
  geom_point() +
  theme_classic() +
  ggrepel::geom_label_repel(data = filter(survival_models, p_val <= sig_level), mapping = aes(label = cluster)) +
  xlab("Distance Between Poor (+) and Favorable (-)") +
  ylab("-log(P value)")
graphics.off()

## Run GSVA on Clusters ----
gene_sets <- list(
  HALLMARK = msigdbr(species = "Homo sapiens", category = "H"),
  GO = msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP"),
  REACTOME = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "REACTOME"),
  KEGG = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG"),
  ONCO = msigdbr(species = "Homo sapiens", category = "C6")
) %>%
  map(~ ListGeneSets(.x, gene_set = "name", gene_name = "symbol"))

gsva_res <- map(gene_sets, ~ RunGSVA(diagnosis, gene_sets = .x, replicates = 3))

stat_res <- map2(gsva_res, names(gsva_res), RunStats, p = 0.05)

pdf(file = here("plots/HALLMARK_GSVA_volcano_plots.pdf"), onefile = TRUE)
stat_res$HALLMARK$plots
graphics.off()

pdf(file = here("plots/GO_GSVA_volcano_plots.pdf"), onefile = TRUE)
stat_res$GO$plots
graphics.off()

pdf(file = here("plots/REACTOME_GSVA_volcano_plots.pdf"), onefile = TRUE)
stat_res$REACTOME$plots
graphics.off()

pdf(file = here("plots/KEGG_GSVA_volcano_plots.pdf"), onefile = TRUE)
stat_res$KEGG$plots
graphics.off()

pdf(file = here("plots/ONCO_GSVA_volcano_plots.pdf"), onefile = TRUE)
stat_res$ONCO$plots
graphics.off()

top_tables <- map(stat_res, ~ .x[["Top Table"]])

top_tables %>%
  reduce(bind_rows) %>%
  list("All Gene Sets" = .) %>%
  append(top_tables) %>%
  openxlsx::write.xlsx(file = here("outs/GSVA_DE_results.xlsx"))

## Run CellPhoneDB ----
if (rerun_cellphone) {
  seurat_tmp <- LoadH5Seurat(file.path(Sys.getenv("AML_DATA"), "05_seurat_annotated.h5Seurat"),
    assays = c("RNA")
  )
  DefaultAssay(seurat_tmp) <- "RNA"
  tmp <- subset(seurat_tmp, timepoint == "Diagnosis" & stemness == "Nonstem")
  tmp <- NormalizeData(tmp)

  Idents(diagnosis) %>%
    as_tibble(rownames = "Cell") %>%
    rename(cell_type = value) %>%
    write_tsv(here("data/cellphonedb_in/metadata.tsv"))

  GetAssayData(tmp, "data") %>%
    DropletUtils::write10xCounts(
      path = here("data/cellphonedb_in/data"),
      overwrite = TRUE
    )

  file.rename(
    here("data/cellphonedb_in/data/genes.tsv"),
    here("data/cellphonedb_in/data/features.tsv")
  )

  rm(tmp, seurat_tmp)

  system(here("lib/cellphonedb.sh"))
}
pdf(file = here("plots/cellphone_plot.pdf"), height = 12, width = 18)
CellPhoneDotPlot(output_path = here("outs/cellphonedb_results"))
graphics.off()

cache_rds(
  expr = {
    cpdb <- ReadCellPhone(output_path = here("outs/cellphonedb_results")) %>%
      pivot_wider(names_from = clusters, values_from = mean) %>%
      select(-pvalue) %>%
      remove_missing(., vars = names(.)[2:length(.)]) %>%
      column_to_rownames(var = "pair") %>%
      CreateSeuratObject()

    cpdb <- ScaleData(cpdb)

    pc <- 100
    cpdb <- RunPCA(cpdb, features = rownames(cpdb), npcs = pc)
    cpdb <- JackStraw(cpdb, dims = pc)
    cpdb <- ScoreJackStraw(cpdb, dims = 1:pc)
    cpdb
  }, hash = file.info(here("outs/cellphonedb_results")),
  file = "cellphone_pca.rds"
)

JackStrawPlot(cpdb, dims = 1:(pc / 2)) + ElbowPlot(cpdb, ndims = (pc / 2)) & NoLegend()

good_pca_vars <- cpdb@reductions$pca@feature.loadings %>%
  as.data.frame() %>%
  abs() %>%
  slice_max(PC_1, n = 50) %>%
  rownames()
pca_vars <- cpdb@reductions$pca@feature.loadings %>%
  as.data.frame() %>%
  .[which(rownames(.) %in% good_pca_vars), ]

iplot <- plot_ly(
  data = pca_vars,
  x = ~PC_1,
  y = ~PC_2,
  text = rownames(pca_vars)
) %>%
  add_segments(x = 0, y = 0, xend = ~PC_1, yend = ~PC_2) %>%
  add_segments(x = 0, xend = 0, y = -0.1, yend = 0.1, color = I("grey"), line = list(dash = "dash")) %>%
  add_segments(y = 0, yend = 0, x = -0.1, xend = 0.1, color = I("grey"), line = list(dash = "dash")) %>%
  layout(
    shapes = list(
      type = "circle",
      xref = "x", x0 = -0.1, x1 = 0.1,
      yref = "y", y0 = -0.1, y1 = 0.1,
      color = I("grey"),
      line = list(dash = "dash")
    ),
    yaxis = list(scaleanchor = "x", zeroline = FALSE, showgrid = FALSE),
    xaxis = list(zeroline = FALSE, showgrid = FALSE),
    showlegend = FALSE,
    title = "CellPhoneDB Feature Loading Plot"
  ) %>%
  add_text(textposition = "top right")

htmlwidgets::saveWidget(as_widget(iplot), here("plots/cellphone_feature_loading.html"))

## DE on Clusters ----
sig_clusters <- survival_models %>%
  filter(chi_sig == "chi_sig") %>%
  pull(cluster)

de_results <- map(sig_clusters, ~ FindMarkers(diagnosis,
  ident.1 = .x,
  test.use = "MAST"
))
names(de_results) <- sig_clusters

de_name <- paste0(sig_clusters, collapse = "+")
de_results[[de_name]] <- FindMarkers(diagnosis,
  ident.1 = sig_clusters,
  test.use = "MAST"
)

de_results[["All Markers"]] <- FindAllMarkers(diagnosis, test.use = "MAST")
openxlsx::write.xlsx(de_results,
  file = here("outs/cluster_DE_results.xlsx"),
  rowNames = TRUE
)

EnhancedVolcano(de_results$`20`, x = "avg_log2FC", y = "adj_p_val")

## Run Granulator ----
### Make References ----
refs <- cache_rds(
  {
    ref <- FindAllMarkers(diagnosis, max.cells.per.ident = 100)
    ref2 <- ref[ref$p_val_adj <= 0.05 & abs(ref$avg_log2FC) >= 0.25, "gene"]
    ref3 <- AverageExpression(diagnosis,
      assays = "SCT",
      features = ref2,
      group.by = "clusters"
    )$SCT

    colnames(ref3) <- colnames(ref3) %>% paste0("cluster", .)

    ciber_gep <- read_tsv(here("cibersort_in/nonstem_clusters_GEP.txt")) %>%
      column_to_rownames(var = "genesymbols") %>%
      as.matrix()

    list(CIBERSORTx = ciber_gep, custom = ref3)
  },
  file = "02-decon_refs.rds",
  hash = list(diagnosis[["clusters"]])
)

sim_plot <- plot_similarity(refs)

### Deconvolute ----

# remove svr if rerunning
# svr is slow and rls seems to be just as effective
methods <- get_decon_methods()[get_decon_methods() %!in% c("svr")] %>% as.vector()
cores <- 12
decon_target <- cache_rds(
  deconvolute(target_data, sigMatrix = refs, use_cores = cores),
  hash = list(target_data, refs),
  file = "02-TARGET_decon.rds"
)
decon_beatAML <- cache_rds(
  deconvolute(beatAML_data, sigMatrix = refs, use_cores = cores),
  hash = list(beatAML_data, refs),
  file = "02-beatAML_decon.rds"
)
decon_tcga <- cache_rds(
  deconvolute(tcga_data, sigMatrix = refs, methods = methods, use_cores = cores),
  hash = list(tcga_data, refs),
  file = "02-TCGA_decon.rds"
)

deconvoluted_samples <- list(
  TARGET = decon_target,
  beatAML = decon_beatAML,
  TCGA = decon_tcga
)
### Make + Print Plots ----
decon_plots <- map(deconvoluted_samples,
  plot_deconvolute,
  scale = TRUE,
  labels = FALSE,
  markers = FALSE
)
cor_res <- map(deconvoluted_samples, correlate)
cor_plots <- map(cor_res,
  plot_correlate,
  method = "heatmap",
  legend = TRUE
)

decon_plots2 <- map(decon_plots, .f = function(ggplot) {
  data <- ggplot$data
  data <- filter(data, celltype == "cluster0")
  ggplot$data <- data
  return(ggplot + facet_wrap(~model))
})

pdf(file = here("plots/decon_plots.pdf"), width = 12, height = 18)
sim_plot
decon_plots
decon_plots2
cor_plots
graphics.off()

### Prepare Data with Clinical Information ----
model_use <- "rls_CIBERSORTx"

target_deconvoluted <- deconvoluted_samples$TARGET$proportions[[model_use]] %>%
  as_tibble(rownames = "patient_USI") %>%
  separate(patient_USI, into = c(NA, NA, "USI", NA, NA), sep = "-") %>%
  left_join(clinical) %>%
  mutate(across(where(is_character), str_to_lower))

target_deconvoluted_subsets <- list(
  FLT3 = filter(target_deconvoluted, `FLT3/ITD positive?` == "yes"),
  NEG = filter(target_deconvoluted, `FLT3/ITD positive?` == "no", `CEBPA mutation` == "no"),
  CEBPA = filter(target_deconvoluted, `CEBPA mutation` == "no")
)

tcga_deconvoluted <- deconvoluted_samples$TCGA$proportions[[model_use]] %>%
  as_tibble(rownames = "case_submitter_id") %>%
  left_join(tcga_ann2) %>%
  filter(flt3_status == "Positive")

beat_aml_deconvoluted <- deconvoluted_samples$beatAML$proportions[[model_use]] %>%
  as_tibble(rownames = "SAMPLE_ID") %>%
  left_join(beat_aml_clinical2) %>%
  filter(FLT3_ITD_CONSENSUS_CALL == "Positive")

## LASSO Model with CIBERSORT ----
sig_clusters <- survival_models %>%
  filter(chi_sig == "chi_sig") %>%
  pull(cluster) %>%
  paste0("cluster", .)

summarise(target_deconvoluted, across(
  .cols = sig_clusters,
  .fns = c(x = mean, med = median)
))

target_deconvoluted_subsets$FLT3 %>%
  mutate(
    status = if_else(`First Event` == "relapse", 1, 0),
    across(
      .cols = starts_with("cluster"),
      .fns = ~ if_else(.x >= mean(.x, na.rm = TRUE), "High", "Low")
    )
  ) %>%
  dplyr::select(
    `Event Free Survival Time in Days`, status,
    all_of(sig_clusters), USI, `FLT3/ITD positive?`
  ) %>%
  remove_missing() %>%
  pivot_longer(cols = starts_with("cluster")) %>%
  group_by(name) %>%
  nest() %>%
  mutate(survival = map(data, ~ survdiff(Surv(`Event Free Survival Time in Days`, status) ~ value, data = .x))) %>%
  mutate(
    chi_sq = map_dbl(survival, ~ .x[["chisq"]]),
    p_val = map_dbl(survival, CalculatePValues.chi)
  ) %>%
  arrange(p_val) %>%
  select(name, chi_sq, p_val) %>%
  write_csv(file = here("outs/target_single_cluster_survival.csv"))

### Train on TARGET (train) ----
target_ciber_split <- split_df(target_deconvoluted_subsets$FLT3,
  y = "Event Free Survival Time in Days",
  ratios = c(0.7, 0.3)
)

lasso_model <- DoLassoModel(target_ciber_split$train,
  `Event Free Survival Time in Days`,
  exclude = c(
    matches("^Year|^Age|^Overall"),
    where(is_character)
  )
)

coef(lasso_model) %>%
  as.data.frame() %>%
  rownames_to_column(var = "predictor") %>%
  rename(coef = s0) %>%
  write_csv(here("outs/TARGET_CIBERSORTx_LASSO_model.csv"))

training <- target_ciber_split$train %>%
  as_tibble() %>%
  mutate(
    status = if_else(`First Event` == "relapse", 1, 0),
    `First Event` = NULL
  ) %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
    score_bin = if_else(score >= median(score), "High", "Low")
  )

target_train <- DoSurvialAnalysis(training,
  `Event Free Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET Training",
  lasso = lasso_model
)

target_train_os <- DoSurvialAnalysis(training,
  `Overall Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET Training",
  lasso = lasso_model
)

### Test on TARGET (test) ----
testing <- target_ciber_split$test %>%
  as_tibble() %>%
  mutate(
    status = if_else(`First Event` == "relapse", 1, 0),
    `First Event` = NULL
  ) %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
    score_bin = if_else(score >= median(score), "High", "Low")
  )

target_test <- DoSurvialAnalysis(testing,
  `Event Free Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET Test",
  lasso = lasso_model
)

target_test_os <- DoSurvialAnalysis(testing,
  `Overall Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET Test",
  lasso = lasso_model
)

### Test on TARGET (FLT3-ITD Negative) ----
neg <- target_deconvoluted_subsets$NEG %>%
  mutate(
    status = if_else(`First Event` == "relapse", 1, 0),
    `First Event` = NULL
  ) %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
    score_bin = if_else(score >= median(score, na.rm = TRUE), "High", "Low")
  )

target_flt3 <- DoSurvialAnalysis(neg,
  `Event Free Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET FLT3-",
  lasso = lasso_model
)

target_flt3_os <- DoSurvialAnalysis(neg,
  `Overall Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET FLT3-",
  lasso = lasso_model
)

### Test on TARGET (CEBPA+) ----
cebpa <- target_deconvoluted_subsets$CEBPA %>%
  mutate(
    status = if_else(`First Event` == "relapse", 1, 0),
    `First Event` = NULL
  ) %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
    score_bin = if_else(score >= median(score, na.rm = TRUE), "High", "Low")
  )

target_cebpa <- DoSurvialAnalysis(cebpa,
  `Event Free Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET CEBPA+",
  lasso = lasso_model
)

target_cebpa_os <- DoSurvialAnalysis(cebpa,
  `Overall Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET CEBPA+",
  lasso = lasso_model
)

### Test on TCGA ----
tcga <- tcga_deconvoluted %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
    score_bin = if_else(score >= median(score), "High", "Low"),
    status = if_else(vital_status == "Dead", 1, 0),
    days_to_death = as.numeric(days_to_death)
  )

tcga_surv <- DoSurvialAnalysis(tcga,
  days_to_death,
  status,
  score,
  group_by = score_bin,
  description = "TCGA (94.1% features present)",
  lasso = lasso_model
)

### Test on BeatAML ----
beat_aml <- beat_aml_deconvoluted %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
    score_bin = if_else(score >= mean(score, na.rm = TRUE), "High", "Low"),
    OS_days = OS_MONTHS * (365.25 / 12)
  )

beat_surv <- DoSurvialAnalysis(beat_aml,
  OS_days,
  status,
  score,
  group_by = score_bin,
  description = "BeatAML",
  lasso = lasso_model
)

### Print Plots ----
pdf(file = here("plots/survival_with_deconv.pdf"))
target_train$plot
target_test$plot
target_train_os$plot
target_test_os$plot
target_flt3$plot
target_flt3_os$plot
target_cebpa$plot
target_cebpa_os$plot
tcga_surv$plot
beat_surv$plot
graphics.off()

## LASSO Model with ONLY CIBERSORTx features ----
target_deconvoluted_subsets$FLT3 %>%
  mutate(
    status = if_else(`First Event` == "relapse", 1, 0),
    across(
      .cols = starts_with("cluster"),
      .fns = ~ if_else(.x >= mean(.x, na.rm = TRUE), "High", "Low")
    )
  ) %>%
  dplyr::select(
    `Event Free Survival Time in Days`, status,
    starts_with("cluster"), USI, `FLT3/ITD positive?`
  ) %>%
  remove_missing() %>%
  pivot_longer(cols = starts_with("cluster")) %>%
  group_by(name) %>%
  nest() %>%
  mutate(survival = map(data, ~ survdiff(Surv(`Event Free Survival Time in Days`, status) ~ value, data = .x))) %>%
  mutate(
    chi_sq = map_dbl(survival, ~ .x[["chisq"]]),
    p_val = map_dbl(survival, CalculatePValues.chi)
  ) %>%
  arrange(p_val) %>%
  select(name, chi_sq, p_val) %>%
  write_csv(file = here("outs/target_single_cluster_survival_CIBERSORTx_only.csv"))

### Train on TARGET (train) ----
lasso_model <- DoLassoModel(target_ciber_split$train,
  `Event Free Survival Time in Days`,
  include = starts_with("cluster")
)

coef(lasso_model) %>%
  as.data.frame() %>%
  rownames_to_column(var = "predictor") %>%
  rename(coef = s0) %>%
  write_csv(here("outs/TARGET_CIBERSORTx_LASSO_model_only_CIBERSORTx.csv"))

training <- target_ciber_split$train %>%
  as_tibble() %>%
  mutate(
    status = if_else(`First Event` == "relapse", 1, 0),
    `First Event` = NULL
  ) %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
    score_bin = if_else(score >= median(score), "High", "Low")
  )

target_train_nf <- DoSurvialAnalysis(training,
  `Event Free Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET Training",
  lasso = lasso_model
)

target_train_os_nf <- DoSurvialAnalysis(training,
  `Overall Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET Training",
  lasso = lasso_model
)

### Test on TARGET (test) ----
testing <- target_ciber_split$test %>%
  as_tibble() %>%
  mutate(
    status = if_else(`First Event` == "relapse", 1, 0),
    `First Event` = NULL
  ) %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
    score_bin = if_else(score >= median(score), "High", "Low")
  )

target_test_nf <- DoSurvialAnalysis(testing,
  `Event Free Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET Test",
  lasso = lasso_model
)

target_test_os_nf <- DoSurvialAnalysis(testing,
  `Overall Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET Training",
  lasso = lasso_model
)

### Test on TARGET (FLT3-ITD Negative) ----
target_ciber_flt3_neg_nf <- target_deconvoluted_subsets$NEG %>%
  mutate(
    status = if_else(`First Event` == "relapse", 1, 0),
    `First Event` = NULL
  ) %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
    score_bin = if_else(score >= median(score, na.rm = TRUE), "High", "Low")
  )

target_flt3_nf <- DoSurvialAnalysis(target_ciber_flt3_neg_nf,
  `Event Free Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET FLT3-",
  lasso = lasso_model
)

target_flt3_os_nf <- DoSurvialAnalysis(target_ciber_flt3_neg_nf,
  `Overall Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET FLT3-",
  lasso = lasso_model
)

### Test on TARGET (CEBPA+) ----
target_ciber_cebpa <- target_deconvoluted_subsets$CEBPA %>%
  mutate(
    status = if_else(`First Event` == "relapse", 1, 0),
    `First Event` = NULL
  ) %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
    score_bin = if_else(score >= median(score, na.rm = TRUE), "High", "Low")
  )

target_cebpa <- DoSurvialAnalysis(target_ciber_cebpa,
  `Event Free Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET CEBPA+",
  lasso = lasso_model
)

target_cebpa_os <- DoSurvialAnalysis(target_ciber_cebpa,
  `Overall Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET CEBPA+",
  lasso = lasso_model
)

### Test on TCGA ----
tcga_ciber <- tcga_deconvoluted %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
    score_bin = if_else(score >= median(score), "High", "Low"),
    status = if_else(vital_status == "Dead", 1, 0),
    days_to_death = as.numeric(days_to_death)
  )

tcga_surv_nf <- DoSurvialAnalysis(tcga_ciber,
  days_to_death,
  status,
  score,
  group_by = score_bin,
  description = "TCGA",
  lasso = lasso_model
)

### Test on BeatAML ----
beat_aml_ciber <- beat_aml_deconvoluted %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
    score_bin = if_else(score >= mean(score, na.rm = TRUE), "High", "Low"),
    OS_days = OS_MONTHS * (365.25 / 12)
  )

beat_surv_nf <- DoSurvialAnalysis(beat_aml_ciber,
  OS_days,
  status,
  score,
  group_by = score_bin,
  description = "BeatAML (67% features present)",
  lasso = lasso_model
)

### Print Plots ----
pdf(file = here("plots/survival_deconv_without_filtering.pdf"))
target_train_nf$plot
target_train_os_nf$plot
target_test_nf$plot
target_test_os_nf$plot
target_flt3_nf$plot
target_flt3_os_nf$plot
target_cebpa$plot
target_cebpa_os$plot
tcga_surv_nf$plot
beat_surv_nf$plot
graphics.off()

## LASSO Model with ONLY CIBERSORTx features ----

### Train on TCGA ----
tcga_ciber <- read_tsv(here("outs/cibersort_results/CIBERSORTx_tcga_Results.txt"))
tcga_ciber <- left_join(tcga_ciber, tcga_ann2, by = c("Mixture" = "case_submitter_id"))

lasso_model_nf <- DoLassoModel(tcga_ciber,
  days_to_death,
  exclude = c(
    matches("^Year|^Age|^Overall|^WBC|%|^lab_|^year_|^age_|cytogen|follow|birth|initial|diagnosis|result"),
    where(is_character),
    where(is_logical),
    where(lubridate::is.instant),
    `P-value`:RMSE,
    Mixture,
    patient_id,
    hydroxyurea_agent_administered_day_count,
    cumulative_agent_total_dose
  )
)

tcga_res <- tcga_ciber %>%
  rename(case_submitter_id = Mixture) %>%
  filter(flt3_status == "Positive") %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(lasso_model_nf)),
    score_bin = if_else(score >= median(score), "High", "Low"),
    status = if_else(vital_status == "Dead", 1, 0),
    days_to_death = as.numeric(days_to_death)
  )

tcga_surv_nf <- DoSurvialAnalysis(tcga_res,
  days_to_death,
  status,
  score,
  group_by = score_bin,
  description = "TCGA",
  lasso = lasso_model_nf
)

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
###   LASSO MODEL SIGNIFICANT   ###
###       FINISH ANALYSIS       ###
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#

## LASSO Model using Features ----
clinical2 <- dplyr::select(
  clinical, `WBC at Diagnosis`, `Bone marrow leukemic blast percentage (%)`,
  `Peripheral blasts (%)`, `Cytogenetic Complexity`,
  `FLT3/ITD allelic ratio`, `Event Free Survival Time in Days`, `FLT3/ITD positive?`, patient_id, `First Event`,
  `Overall Survival Time in Days`
)

gep <- read_tsv(here("cibersort_in/nonstem_GEP.txt"))

genes <- gep %>%
  dplyr::select(cluster18, cluster23, genesymbols) %>%
  filter(cluster18 != 1 | cluster23 != 1) %>%
  filter(cluster18 == 1 | cluster23 == 1) %>%
  pull(genesymbols)

target <- read_tsv(file = here("cibersort_in/target_data.txt"))

target_before_split <- target %>%
  transposeDF() %>%
  dplyr::select(samples, any_of(genes)) %>%
  mutate(across(.cols = 2:ncol(.), .fns = as.numeric)) %>%
  separate(col = "samples", into = c(NA, NA, "patient_id", NA, NA)) %>%
  inner_join(clinical2) %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "yes") %>%
  dplyr::select(-`FLT3/ITD positive?`, -patient_id) %>%
  mutate(across(.cols = !`First Event`, .fns = as.numeric)) %>%
  remove_missing()

target_split <- split_df(target_before_split, y = "Event Free Survival Time in Days", ratios = c(0.7, 0.3))

### Train on TARGET (train) ----
target_model <- DoLassoModel(target_split$train,
  outcome = `Event Free Survival Time in Days`,
  exclude = c(`First Event`, `Overall Survival Time in Days`)
)

coef(target_model)

target_train_score <- target_split$train %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(target_model)),
    score_bin = if_else(score >= median(score), "High", "Low"),
    status = if_else(`First Event` == "relapse", 1, 0)
  )

target_train_surv <- DoSurvialAnalysis(target_train_score,
  time = `Event Free Survival Time in Days`,
  status = status,
  predictor = score,
  group_by = score_bin,
  description = "TARGET Training Data",
  lasso = target_model
)

target_train_surv_os <- DoSurvialAnalysis(target_train_score,
  time = `Overall Survival Time in Days`,
  status = status,
  predictor = score,
  group_by = score_bin,
  description = "TARGET Training Data",
  lasso = target_model
)

### Test with TARGET (test) ----
target_test_score <- target_split$test %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(target_model)),
    score_bin = if_else(score >= median(score), "High", "Low"),
    status = if_else(`First Event` == "relapse", 1, 0)
  )

target_test_surv <- DoSurvialAnalysis(target_test_score,
  time = `Event Free Survival Time in Days`,
  status = status,
  predictor = score,
  group_by = score_bin,
  description = "TARGET Test Data",
  lasso = target_model
)

target_test_surv_os <- DoSurvialAnalysis(target_test_score,
  time = `Overall Survival Time in Days`,
  status = status,
  predictor = score,
  group_by = score_bin,
  description = "TARGET Test Data",
  lasso = target_model
)

### Test on TARGET (FLT3-ITD Negative) ----
target_ciber_flt3_neg <- target_ciber %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`FLT3/ITD positive?` == "no" & `CEBPA mutation` == "no")

target_ciber_flt3_neg_feats <- target_ciber_flt3_neg %>%
  mutate(
    status = if_else(`First Event` == "relapse", 1, 0),
    `First Event` = NULL
  ) %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
    score_bin = if_else(score >= median(score, na.rm = TRUE), "High", "Low")
  )

target_flt3_feats <- DoSurvialAnalysis(target_ciber_flt3_neg_feats,
  `Event Free Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET FLT3-",
  lasso = lasso_model
)

target_flt3_os_feats <- DoSurvialAnalysis(target_ciber_flt3_neg_feats,
  `Overall Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET FLT3-",
  lasso = lasso_model
)

### Test on TARGET (CEBPA+) ----
target_ciber_cebpa <- target_ciber %>%
  mutate(across(where(is_character), str_to_lower)) %>%
  filter(`CEBPA mutation` == "yes")

target_ciber_cebpa <- target_ciber_cebpa %>%
  mutate(
    status = if_else(`First Event` == "relapse", 1, 0),
    `First Event` = NULL
  ) %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
    score_bin = if_else(score >= median(score, na.rm = TRUE), "High", "Low")
  )

target_cebpa <- DoSurvialAnalysis(target_ciber_cebpa,
  `Event Free Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET CEBPA+",
  lasso = lasso_model
)

target_cebpa_os <- DoSurvialAnalysis(target_ciber_cebpa,
  `Overall Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET CEBPA+",
  lasso = lasso_model
)

### Test with TCGA on TARGET model ----
tcga_data <- read_tsv(here("cibersort_in/tcga_cibersort.txt"))
tcga_data <- transposeDF(tcga_data, rowname_col = "case_submitter_id")

annotated_tcga <- read_tsv(here("data/tcga/clinical.cart.2021-04-22/clinical.tsv"), na = "'--") %>%
  right_join(tcga_data) %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(target_model)),
    score_bin = if_else(score >= median(score), "High", "Low"),
    status = if_else(vital_status == "Dead", 1, 0)
  )

tcga_surv <- DoSurvialAnalysis(annotated_tcga,
  time = days_to_death,
  status = status,
  predictor = score,
  group_by = score_bin,
  description = "TCGA Data (4.8% features present)",
  lasso = target_model
)

### Test with BeatAML on TARGET model ----
beat_aml <- read_tsv(here("data/aml_ohsu_2018/data_RNA_Seq_expression_cpm.txt"))

beat_aml_data <- beat_aml %>%
  dplyr::select(-Entrez_Gene_Id) %>%
  transposeDF(rowname_col = SAMPLE_ID) %>%
  mutate(across(.cols = !starts_with("SAMPLE_ID"), .fns = as.numeric)) %>%
  inner_join(beat_aml_clinical2) %>%
  filter(FLT3_ITD_CONSENSUS_CALL == "Positive") %>%
  mutate(across(.cols = !any_of(c("FLT3_ITD_CONSENSUS_CALL", "SAMPLE_ID")), .fns = as.numeric),
    OS_DAYS = OS_MONTHS * 365.25 / 12,
    score = UseLASSOModelCoefs(cur_data(), coef(target_model)),
    score_bin = if_else(score >= mean(score), "High", "Low")
  )

beataml_surv <- DoSurvialAnalysis(beat_aml_data,
  time = OS_DAYS,
  status = status,
  predictor = score,
  group_by = score_bin,
  description = "Beat-AML Data (85.7% features present)",
  lasso = target_model
)

### Print Plots ----
pdf(here("plots/features_EFS_survival.pdf"))
target_train_surv$plot
target_train_surv_os$plot
target_test_surv$plot
target_test_surv_os$plot
target_test_arm$plot
target_cebpa$plot
target_cebpa_os$plot
tcga_surv$plot
beataml_surv$plot
graphics.off()

## pLSC6 score ----
pLSC6 <- matrix(
  data = c(0, 0.189, 0.054, 0.0171, 0.141, 0.109, 0.0516),
  ncol = 1,
  dimnames = list(
    c("(Intercept)", "DNMT3B", "GPR56", "CD34", "SOCS2", "SPINK2", "FAM30A"),
    c("s0")
  )
) %>% as.sparse()

### Test with TARGET ----
target_clinical <- target_data %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  transposeDF(colname_col = "rowname") %>%
  separate(col = samples, into = c(NA, NA, "USI", NA, NA)) %>%
  inner_join(clinical)

target_clinical_subset <- list(
  FLT = filter(
    target_clinical,
    `FLT3/ITD positive?` == "Yes"
  ),
  CEBPA = filter(
    target_clinical,
    `CEBPA mutation` == "Yes"
  ),
  NEG = filter(
    target_clinical,
    `CEBPA mutation` == "No",
    `FLT3/ITD positive?` == "No"
  )
)

target_test_score <- target_clinical_subset$FLT %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), pLSC6),
    score_bin = if_else(score >= median(score), "High", "Low"),
    status = if_else(`First Event` == "relapse", 1, 0)
  )

target_test_surv <- DoSurvialAnalysis(target_test_score,
  time = `Event Free Survival Time in Days`,
  status = status,
  predictor = score,
  group_by = score_bin,
  description = "TARGET Test Data"
)

target_test_surv_os <- DoSurvialAnalysis(target_test_score,
  time = `Overall Survival Time in Days`,
  status = status,
  predictor = score,
  group_by = score_bin,
  description = "TARGET Test Data"
)

### Test on TARGET (FLT3-ITD Negative) ----
target_ciber_flt3_neg_feats <- target_clinical_subset$NEG %>%
  mutate(
    status = if_else(`First Event` == "relapse", 1, 0),
    `First Event` = NULL
  ) %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), pLSC6),
    score_bin = if_else(score >= median(score, na.rm = TRUE), "High", "Low")
  )

target_flt3_feats <- DoSurvialAnalysis(target_ciber_flt3_neg_feats,
  `Event Free Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET FLT3-"
)

target_flt3_os_feats <- DoSurvialAnalysis(target_ciber_flt3_neg_feats,
  `Overall Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET FLT3-"
)

### Test on TARGET (CEBPA+) ----
target_ciber_cebpa <- target_clinical_subset$CEBPA %>%
  mutate(
    status = if_else(`First Event` == "relapse", 1, 0),
    `First Event` = NULL
  ) %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), pLSC6),
    score_bin = if_else(score >= median(score, na.rm = TRUE), "High", "Low")
  )

target_cebpa <- DoSurvialAnalysis(target_ciber_cebpa,
  `Event Free Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET CEBPA+"
)

target_cebpa_os <- DoSurvialAnalysis(target_ciber_cebpa,
  `Overall Survival Time in Days`,
  status,
  score,
  group_by = score_bin,
  description = "TARGET CEBPA+"
)

### Test on TCGA ----
tcga_ciber <- tcga_data %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  transposeDF(colname_col = "rowname", rowname_col = "case_submitter_id") %>%
  inner_join(tcga_ann2) %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
    score_bin = if_else(score >= median(score), "High", "Low"),
    status = if_else(vital_status == "Dead", 1, 0),
    days_to_death = as.numeric(days_to_death)
  )

tcga_surv_nf <- DoSurvialAnalysis(tcga_ciber,
  days_to_death,
  status,
  score,
  group_by = score_bin,
  description = "TCGA",
  lasso = lasso_model
)

### Test on BeatAML ----
beat_aml_ciber <- beatAML_data %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  transposeDF(rowname_col = "SAMPLE_ID", colname_col = "rowname") %>%
  inner_join(beat_aml_clinical2) %>%
  mutate(
    score = UseLASSOModelCoefs(cur_data(), coef(lasso_model)),
    score_bin = if_else(score >= mean(score, na.rm = TRUE), "High", "Low"),
    OS_days = OS_MONTHS * (365.25 / 12)
  )

beat_surv_nf <- DoSurvialAnalysis(beat_aml_ciber,
  OS_days,
  status,
  score,
  group_by = score_bin,
  description = "BeatAML (67% features present)",
  lasso = lasso_model
)

### Print Plots ----
pdf(here("plots/pLSC6_survival.pdf"))
target_test_surv$plot
target_test_surv_os$plot
target_flt3_feats$plot
target_flt3_os_feats$plot
target_cebpa$plot
target_cebpa_os$plot
graphics.off()
