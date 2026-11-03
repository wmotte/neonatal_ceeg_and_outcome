#!/usr/bin/env Rscript
#
# lib_ml.R
# Leakage-free prediction pipeline.
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Every data-dependent step is learned on the training part of a split only and
# then applied to the held-out part:
#   1. median (continuous) / mode (binary) imputation
#   2. optional PMA-normalisation of quantitative cEEG features (residuals of a
#      training-fold linear regression on PMA)
#   3. hyperparameter tuning (inner stratified k-fold, log-loss)
#   4. classification threshold (= outcome prevalence in the training part)
# Outer validation: stratified 5-fold cross-validation, repeated 10 times.
# Algorithms: LightGBM (primary, as in the original submission) and
# ridge-penalised logistic regression (comparator).
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lightgbm)
  library(glmnet)
})

# ---- Helpers -----------------------------------------------------------------
logit  <- function(p) log(p / (1 - p))
clip01 <- function(p, eps = 1e-4) pmin(pmax(p, eps), 1 - eps)

# Fast Mann-Whitney AUC.
auc_fast <- function(y, p) {
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Stratified fold labels (1..k), random within each outcome class.
strat_folds <- function(y, k) {
  f <- integer(length(y))
  for (cl in unique(y)) {
    idx <- which(y == cl)
    f[idx] <- sample(rep_len(seq_len(k), length(idx)))
  }
  f
}

# ---- Encoding ----------------------------------------------------------------
# Factors -> 0/1 dummies (first level = reference). Returns a numeric data
# frame plus a column -> original-variable map (to aggregate SHAP values).
encode <- function(df, vars) {
  out <- list(); map <- character()
  for (v in vars) {
    x <- df[[v]]
    if (is.factor(x)) {
      for (lv in levels(x)[-1]) {
        nm <- paste0(v, "__", make.names(lv))
        out[[nm]] <- as.numeric(x == lv)
        map[nm] <- v
      }
    } else {
      out[[v]] <- as.numeric(x); map[v] <- v
    }
  }
  list(X = as_tibble(out), map = map)
}

# ---- Preprocessing learned on training data only -----------------------------
prep_learn <- function(Xtr, pma_tr = NULL, normalise_cols = character()) {
  binary <- map_lgl(Xtr, ~ all(na.omit(.x) %in% c(0, 1)))
  fill <- map2_dbl(Xtr, binary, function(x, b) {
    if (b) as.numeric(mean(x, na.rm = TRUE) >= 0.5) else median(x, na.rm = TRUE)
  })
  norm <- list()
  if (length(normalise_cols) && !is.null(pma_tr)) {
    Ximp <- impute_apply(Xtr, fill)
    for (cc in intersect(normalise_cols, names(Xtr))) {
      norm[[cc]] <- coef(lm(Ximp[[cc]] ~ pma_tr))
    }
  }
  list(fill = fill, norm = norm)
}
impute_apply <- function(X, fill) {
  for (cc in names(fill)) X[[cc]][is.na(X[[cc]])] <- fill[[cc]]
  X
}
prep_apply <- function(rec, X, pma = NULL) {
  X <- impute_apply(X, rec$fill)
  for (cc in names(rec$norm)) {
    b <- rec$norm[[cc]]
    X[[cc]] <- X[[cc]] - (b[1] + b[2] * pma)
  }
  as.matrix(X)
}

# ---- LightGBM with inner-CV tuning -------------------------------------------
lgb_params <- function(g, seed) {
  list(objective = "binary", learning_rate = g$learning_rate,
       max_depth = g$max_depth, num_leaves = g$num_leaves,
       min_data_in_leaf = g$min_data_in_leaf, min_sum_hessian_in_leaf = 1e-3,
       feature_pre_filter = FALSE, force_col_wise = TRUE, num_threads = 1L,
       deterministic = TRUE, seed = seed, verbose = -1L)
}
lgb_fit_one <- function(X, y, g, nrounds, seed) {
  lgb.train(params = lgb_params(g, seed),
            data = lgb.Dataset(X, label = y, params = list(feature_pre_filter = FALSE)),
            nrounds = nrounds, verbose = -1L)
}
fit_lgb <- function(X, y, grid = CFG$lgb_grid, inner = CFG$inner_folds, seed = 1L) {
  # Tune max_depth x min_data_in_leaf by fitting the largest nrounds once per
  # inner fold and scoring every candidate nrounds from the same booster.
  base <- distinct(grid, max_depth, num_leaves, min_data_in_leaf, learning_rate)
  nr   <- sort(unique(grid$nrounds))
  fid  <- strat_folds(y, inner)
  loss <- matrix(0, nrow(base), length(nr))
  for (b in seq_len(nrow(base))) {
    for (k in seq_len(inner)) {
      tr <- fid != k
      m <- lgb_fit_one(X[tr, , drop = FALSE], y[tr], base[b, ], max(nr), seed)
      for (j in seq_along(nr)) {
        p <- clip01(predict(m, X[!tr, , drop = FALSE], num_iteration = nr[j]))
        yy <- y[!tr]
        loss[b, j] <- loss[b, j] - sum(yy * log(p) + (1 - yy) * log(1 - p))
      }
    }
  }
  best <- which(loss == min(loss), arr.ind = TRUE)[1, ]
  g <- base[best[1], ]; n_best <- nr[best[2]]
  list(algo = "lgb", model = lgb_fit_one(X, y, g, n_best, seed),
       tuned = mutate(g, nrounds = n_best))
}

# ---- Ridge logistic regression -----------------------------------------------
fit_ridge <- function(X, y, inner = CFG$inner_folds, seed = 1L) {
  keep <- apply(X, 2, function(z) sd(z) > 0)
  fid <- strat_folds(y, inner)
  cv <- cv.glmnet(X[, keep, drop = FALSE], y, family = "binomial", alpha = 0,
                  foldid = fid, type.measure = "deviance", standardize = TRUE)
  list(algo = "ridge", model = cv, keep = keep, tuned = tibble(lambda = cv$lambda.min))
}

predict_fit <- function(fit, X) {
  if (fit$algo == "lgb") return(predict(fit$model, X))
  as.numeric(predict(fit$model, X[, fit$keep, drop = FALSE], s = "lambda.min", type = "response"))
}

# TreeSHAP contributions (log-odds scale), aggregated to original variables.
shap_fit <- function(fit, X, map) {
  if (fit$algo != "lgb") return(NULL)
  cmat <- predict(fit$model, X, type = "contrib")
  cmat <- cmat[, seq_len(ncol(X)), drop = FALSE]
  colnames(cmat) <- colnames(X)
  vars <- unique(map[colnames(X)])
  sapply(vars, function(v) rowSums(cmat[, map[colnames(X)] == v, drop = FALSE]))
}

# ---- Model specification -----------------------------------------------------
# spec: list(features = <vars>, algo = "lgb"|"ridge", pma_normalise = FALSE)
# The PMA-normalisation only touches the quantitative cEEG features.
build_X <- function(df, spec) encode(df, spec$features)

# ---- Repeated stratified CV ---------------------------------------------------
# Returns out-of-fold predictions (+ fold thresholds) and per-fold SHAP summaries.
run_cv <- function(df, outcome, spec, repeats = CFG$cv_repeats, folds = CFG$cv_folds,
                   seed = CFG$seed, keep_shap = TRUE) {
  y   <- df[[outcome]]
  enc <- build_X(df, spec)
  X0  <- enc$X
  norm_cols <- if (isTRUE(spec$pma_normalise)) intersect(names(CFG$eeg_quant), names(X0)) else character()
  oof <- vector("list", repeats * folds); shp <- vector("list", repeats * folds); i <- 0
  for (r in seq_len(repeats)) {
    set.seed(seed + r)
    fid <- strat_folds(y, folds)
    for (k in seq_len(folds)) {
      i <- i + 1
      tr <- fid != k; te <- !tr
      rec <- prep_learn(X0[tr, ], df$pma[tr], norm_cols)
      Xtr <- prep_apply(rec, X0[tr, ], df$pma[tr])
      Xte <- prep_apply(rec, X0[te, ], df$pma[te])
      fit <- if (spec$algo == "lgb") fit_lgb(Xtr, y[tr], seed = seed + i) else fit_ridge(Xtr, y[tr], seed = seed + i)
      p   <- predict_fit(fit, Xte)
      thr <- mean(y[tr])
      oof[[i]] <- tibble(rep = r, fold = k, row = which(te), record = df$record[te],
                         y = y[te], prob = p, thr = thr)
      if (keep_shap && spec$algo == "lgb") {
        s <- shap_fit(fit, Xte, enc$map)
        shp[[i]] <- tibble(rep = r, fold = k, variable = colnames(s),
                           mean_abs_shap = colMeans(abs(s)))
      }
    }
  }
  list(oof = bind_rows(oof), shap = bind_rows(shp))
}

# ---- Performance metrics -----------------------------------------------------
metrics_one <- function(y, p, thr) {
  cls  <- as.integer(p >= thr)
  tp <- sum(cls == 1 & y == 1); tn <- sum(cls == 0 & y == 0)
  fp <- sum(cls == 1 & y == 0); fn <- sum(cls == 0 & y == 1)
  sens <- tp / (tp + fn); spec <- tn / (tn + fp)
  pc <- clip01(p); lp <- logit(pc)
  cal_slope <- tryCatch(unname(coef(glm(y ~ lp, family = binomial))[2]), error = function(e) NA_real_)
  cal_int   <- tryCatch(unname(coef(glm(y ~ 1 + offset(lp), family = binomial))[1]), error = function(e) NA_real_)
  c(auc = auc_fast(y, p), bacc = (sens + spec) / 2, sens = sens, spec = spec,
    ppv = if ((tp + fp) > 0) tp / (tp + fp) else NA_real_,
    npv = if ((tn + fn) > 0) tn / (tn + fn) else NA_real_,
    brier = mean((p - y)^2), cal_intercept = cal_int, cal_slope = cal_slope)
}

# Point estimate = mean over CV repeats. 95% CI = percentile bootstrap over
# subjects (each bootstrap sample re-evaluates every repeat, then averages).
summarise_cv <- function(oof, n_boot = CFG$n_boot, seed = CFG$seed, ci_metrics = c("auc", "bacc", "sens", "spec", "ppv", "brier")) {
  reps <- split(oof, oof$rep)
  ids  <- sort(unique(oof$row))
  point <- reps |> map(~ metrics_one(.x$y, .x$prob, .x$thr)) |> reduce(`+`) / length(reps)
  set.seed(seed)
  # pre-index each repeat by subject row for fast resampling
  reps_m <- map(reps, ~ arrange(.x, row))
  boot <- replicate(n_boot, {
    b <- sample(seq_along(ids), replace = TRUE)
    m <- map(reps_m, function(rr) {
      yy <- rr$y[b]; pp <- rr$prob[b]; tt <- rr$thr[b]
      cls <- as.integer(pp >= tt)
      sens <- sum(cls == 1 & yy == 1) / sum(yy == 1); spec <- sum(cls == 0 & yy == 0) / sum(yy == 0)
      c(auc = auc_fast(yy, pp), bacc = (sens + spec) / 2, sens = sens, spec = spec,
        ppv = sum(cls == 1 & yy == 1) / max(1, sum(cls == 1)), brier = mean((pp - yy)^2))
    })
    reduce(m, `+`) / length(m)
  })
  ci <- apply(boot, 1, quantile, c(.025, .975), na.rm = TRUE)
  tibble(metric = names(point), estimate = unname(point)) |>
    left_join(tibble(metric = colnames(ci), lo = ci[1, ], hi = ci[2, ]), by = "metric")
}

# Paired bootstrap for the AUC difference between two models evaluated on the
# same outer splits (same seed -> identical fold assignments).
delta_auc <- function(oof_a, oof_b, n_boot = CFG$n_boot, seed = CFG$seed) {
  a <- split(arrange(oof_a, rep, row), oof_a$rep |> sort())
  b <- split(arrange(oof_b, rep, row), oof_b$rep |> sort())
  est <- mean(map2_dbl(a, b, ~ auc_fast(.x$y, .x$prob) - auc_fast(.y$y, .y$prob)))
  n <- nrow(a[[1]]); set.seed(seed)
  bt <- replicate(n_boot, {
    s <- sample.int(n, replace = TRUE)
    mean(map2_dbl(a, b, ~ auc_fast(.x$y[s], .x$prob[s]) - auc_fast(.y$y[s], .y$prob[s])))
  })
  tibble(delta_auc = est, lo = quantile(bt, .025), hi = quantile(bt, .975),
         p_boot = 2 * min(mean(bt <= 0), mean(bt >= 0)))
}

# ---- Permutation tests (complete pipeline re-run for every permutation) ------
# type = "label": permute the outcome (null: no association at all).
# type = "eeg"  : permute the cEEG block across subjects, keep clinical data
#                 and outcome intact (null: cEEG adds nothing beyond clinical).
perm_test <- function(df, outcome, spec, type = c("label", "eeg"), n_perm = CFG$n_perm,
                      seed = CFG$seed) {
  type <- match.arg(type)
  eeg_vars <- intersect(names(CFG$eeg_all), spec$features)
  furrr::future_map_dbl(seq_len(n_perm), function(b) {
    set.seed(seed + 100000L + b)
    dp <- df
    perm <- sample.int(nrow(df))
    if (type == "label") dp[[outcome]] <- dp[[outcome]][perm]
    else dp[eeg_vars] <- dp[perm, eeg_vars]
    res <- run_cv(dp, outcome, spec, repeats = 1L, seed = seed + 200000L + b, keep_shap = FALSE)
    auc_fast(res$oof$y, res$oof$prob)
  }, .options = furrr::furrr_options(seed = TRUE, packages = c("tidyverse", "lightgbm", "glmnet")))
}
