suppressPackageStartupMessages({
  library(here)
  library(log4r)
  library(survival)
  library(survminer)
  library(scorecard)
  library(glmnet)
  library(tidyverse)
  library(glue)
})

UseLASSOModelCoefs <- function(x, coef) {
  mod2 <- coef[-1, , drop = FALSE] #%>% as.matrix() %>% as.data.frame()

  w_vec <- list(0)
  mod2 <- mod2[which(mod2 != 0), ,drop = FALSE]

  missing_features <- mod2[which(rownames(mod2) %!in% colnames(x)), , drop = FALSE] %>% rownames()
  n_features <- nrow(mod2)
  n_present <- n_features - length(missing_features)

  if (n_present == 0) rlang::abort(message = "There are no features matching the lasso model present!")

  if (length(missing_features) >= 1)  {
    # build message
    mfeat <- "The following features are missing in the data: "
    feats <- paste(missing_features, collapse = ", ")
    pfeat <- paste0("Only ",
                    round(n_present / n_features * 100, digits = 1),
                    "% features present.")
    rlang::warn(message = c(mfeat, feats, pfeat))
  }
  for (i in seq_along(rownames(mod2))) {
    col <- rownames(mod2)[i]
    vec <- x[[col]]
    w <- mod2[i, 1]
    w_vec[[i]] <- vec * w
  }
  w_vec <- purrr::discard(w_vec, is_empty)
  w_len <- length(w_vec)

  if (w_len == 1) {
    res <- w_vec[[1]]
  } else {
    res <- w_vec %>% purrr::reduce(rbind) %>% colSums()
  }
  return(res + coef[1, ])
}

`%!in%` <- function(lhs, rhs) {
  match(lhs, rhs, nomatch = 0) == 0
}

# stem Model: outs/flt3_cebpa_stem/EFS_with_clinical_feats/glm_coef_alpha_0.25.tsv
# nonstem model: outs/flt3_cebpa/EFS_with_clinical/glm_coef_alpha_0.25.tsv
lasso_model <- read_tsv("outs/flt3_cebpa/EFS_with_clinical/glm_coef_alpha_0.25.tsv")
names(lasso_model[[2]]) <- lasso_model[[1]]
lasso_model <- lasso_model[[2]]
#lasso_model <- readRDS("outs/all_cells/EFS/lasso_model.rds")
data <- readRDS("outs/all_cells/cache/clinical_deconvoluted.rds")

lasso_model_red <- lasso_model %>%
  as.matrix() %>%
  as.data.frame()

lasso_model_red <- as.data.frame(lasso_model_red)

training <- data[["TRAIN"]]
training[["score"]] <- UseLASSOModelCoefs(training, lasso_model_red)
rmse <- sqrt(sum((training[["score"]] - training[["Event Free Survival Time in Days"]])^2) / nrow(training))

# dicotiomize score -> use this for parameters to assess
flt3 <- data[["TARGET:FLT3"]]
flt3[["score"]] <- UseLASSOModelCoefs(flt3, lasso_model_red)
# high surv = 0
flt3 %>%
  mutate(
    score_bin = if_else(score > median(score), 0, 1),
    efs_bin = if_else(`Event Free Survival Time in Days` > median(`Event Free Survival Time in Days`), 0, 1),
    false_neg = if_else(efs_bin == 1 & score_bin == 0, 1, 0),
    false_pos = if_else(efs_bin == 0 & score_bin == 1, 1, 0),
    true_pos = if_else(efs_bin == 1 & score_bin == 1, 1, 0),
    true_neg = if_else(efs_bin == 0 & score_bin == 0, 1, 0)
  ) %>%
  summarise(
    accuracy = (sum(true_pos) + sum(true_neg)) / nrow(flt3),
    sensitivity = sum(score_bin) / (sum(efs_bin) + sum(false_neg)),
    specificity = (nrow(flt3) - sum(efs_bin)) / (sum(false_pos) + nrow(flt3) - sum(efs_bin))
  )


# dicotiomize score using cutp

