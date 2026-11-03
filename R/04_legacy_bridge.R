#!/usr/bin/env Rscript
#
# 04_legacy_bridge.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Bridging analysis between the initial pipeline (X. Wang, Python notebook in
# data/legacy/) and the current pipeline.
#
# What the notebook did (verified line by line):
#   * median imputation of the quantitative features on the FULL cohort;
#   * a "propensity score" computed in SPSS on the FULL cohort. We verified it
#     is the fitted probability of the OUTCOME from a logistic regression on GA,
#     neonatal seizures (4 categories), brain injury and medication
#     (max abs difference < 1e-9), i.e. a prognostic (disease-risk) score;
#   * OLS residualisation of all 14 features on that score + standardisation,
#     on the FULL cohort;
#   * LightGBM (dart, max_depth = 1, 200 trees, lr 0.1, min_child_samples 20;
#     hyperparameters read from the pickled models), tuned outside the CV loop;
#   * leave-one-out CV with the default 0.5 threshold; no permutation code and
#     no per-fold SHAP code in the notebook.
#
# Here we (L1) re-create that pipeline in R and (L2) re-run it with every step
# moved inside the LOO training folds. We also show why the binary asymmetry
# feature became the top SHAP feature after residualisation.
# ==============================================================================

if (!exists("CFG")) source(file.path("R", "00_config.R"))
if (!exists("run_cv")) source(file.path(CFG$root, "R", "lib_ml.R"))
dat <- readRDS(file.path(CFG$dir_derived, "analysis_data.rds"))

feat  <- names(CFG$eeg_all)
legacy_par <- list(objective = "binary", boosting = "dart", max_depth = 1L, num_leaves = 31L,
                   learning_rate = 0.1, min_data_in_leaf = 20L, min_sum_hessian_in_leaf = 1e-3,
                   feature_pre_filter = FALSE, force_col_wise = TRUE, num_threads = 1L,
                   deterministic = TRUE, seed = CFG$seed, verbose = -1L)
fit_legacy <- function(X, y) lgb.train(legacy_par, lgb.Dataset(X, label = y), nrounds = 200L, verbose = -1L)

ps_formula <- function(outcome) reformulate(c("ga", "seizure_cat", "brain_injury", "medication"), outcome)

residualise <- function(Xtr, ps_tr, Xte = NULL, ps_te = NULL) {
  B <- sapply(seq_len(ncol(Xtr)), function(j) coef(lm(Xtr[, j] ~ ps_tr)))
  rtr <- Xtr - cbind(1, ps_tr) %*% B
  mu <- colMeans(rtr); s <- apply(rtr, 2, sd); s[s == 0] <- 1
  out <- list(tr = sweep(sweep(rtr, 2, mu), 2, s, "/"))
  if (!is.null(Xte)) out$te <- sweep(sweep(Xte - cbind(1, ps_te) %*% B, 2, mu), 2, s, "/")
  out
}

legacy_metrics <- function(y, p, thr = 0.5) {
  m <- metrics_one(y, p, thr)
  m[c("auc", "bacc", "sens", "spec", "ppv")]
}

res <- list(); shap_leg <- list(); sym_demo <- list()
for (o in names(CFG$outcomes)) {
  d <- filter(dat, !is.na(.data[[o]]))
  y <- d[[o]]
  X <- as.matrix(d[feat])

  # ---- L1: legacy (full-cohort preprocessing) ----
  rec_full <- prep_learn(as_tibble(X))
  Xi <- prep_apply(rec_full, as_tibble(X))
  ps_full <- predict(glm(ps_formula(o), data = d, family = binomial), type = "response")
  Xr <- residualise(Xi, ps_full)$tr
  p1 <- map_dbl(seq_along(y), function(i) predict(fit_legacy(Xr[-i, ], y[-i]), Xr[i, , drop = FALSE]))

  # ---- L2: identical model, every step inside the LOO training fold ----
  p2 <- map_dbl(seq_along(y), function(i) {
    tr <- setdiff(seq_along(y), i)
    rec <- prep_learn(as_tibble(X[tr, ]))
    Xtr <- prep_apply(rec, as_tibble(X[tr, ])); Xte <- prep_apply(rec, as_tibble(X[i, , drop = FALSE]))
    g <- glm(ps_formula(o), data = d[tr, ], family = binomial)
    pstr <- predict(g, type = "response"); pste <- predict(g, newdata = d[i, ], type = "response")
    rr <- residualise(Xtr, pstr, Xte, pste)
    predict(fit_legacy(rr$tr, y[tr]), rr$te)
  })

  res[[o]] <- bind_rows(
    c(outcome = o, pipeline = "L1 legacy: preprocessing on full cohort", legacy_metrics(y, p1)),
    c(outcome = o, pipeline = "L2 legacy model, preprocessing inside folds", legacy_metrics(y, p2)))

  # ---- SHAP of the legacy model fitted on the full cohort ----
  m_full <- fit_legacy(Xr, y)
  cm <- predict(m_full, Xr, type = "contrib")[, seq_along(feat)]
  shap_leg[[o]] <- tibble(outcome = o, variable = feat, mean_abs_shap = colMeans(abs(cm)))

  # ---- Why asymmetry: after residualisation, the binary asymmetry column is a
  #      linear function of the prognostic score in the 93% symmetric infants ----
  sym_demo[[o]] <- tibble(outcome = o, record = d$record, asymmetry = Xi[, "asymmetry"],
                          ps = ps_full, asym_resid = Xr[, "asymmetry"], y = y)
}

legacy <- bind_rows(res) |> mutate(across(c(auc, bacc, sens, spec, ppv), as.numeric))
write_tsv(legacy, file.path(CFG$dir_tab, "tab_legacy_bridge.tsv"))
shap_leg <- bind_rows(shap_leg) |> group_by(outcome) |>
  mutate(rank = rank(-mean_abs_shap)) |> ungroup()
write_tsv(shap_leg, file.path(CFG$dir_tab, "tab_legacy_shap.tsv"))
sym_demo <- bind_rows(sym_demo)
sym_cor <- sym_demo |> group_by(outcome) |>
  summarise(cor_resid_ps_symmetric = cor(asym_resid[asymmetry == 0], ps[asymmetry == 0]),
            cor_resid_ps_all = cor(asym_resid, ps),
            auc_asym_resid = auc_fast(y, -asym_resid), auc_asym_raw = auc_fast(y, asymmetry),
            auc_ps = auc_fast(y, ps), .groups = "drop")
write_tsv(sym_cor, file.path(CFG$dir_tab, "tab_legacy_asymmetry_artefact.tsv"))
saveRDS(list(legacy = legacy, shap = shap_leg, sym = sym_demo), file.path(CFG$dir_derived, "legacy.rds"))
message("[04] legacy bridge done")
