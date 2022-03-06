UseLASSOModelCoefs <- function(x, coef) {
  mod2 <- coef[-1, ] %>% as.data.frame()
  w_vec <- list(0)
  mod2 <- mod2[which(mod2[1] != 0), , drop = FALSE]
  missing_features <- mod2[which(rownames(mod2) %!in% colnames(x)), ,
    drop = FALSE
  ] %>% rownames()
  n_features <- nrow(mod2)
  n_present <- n_features - length(missing_features)
  if (n_present == 0) {
    rlang::abort(message = "There are no features matching the lasso model present!")
  }
  if (length(missing_features) >= 1) {
    mfeat <- "The following features are missing in the data: "
    feats <- paste(missing_features, collapse = ", ")
    pfeat <- paste0("Only ", round(n_present / n_features *
      100, digits = 1), "% features present.")
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
    res <- w_vec %>%
      purrr::reduce(rbind) %>%
      colSums()
  }
  return(res + coef[1, ])
}
