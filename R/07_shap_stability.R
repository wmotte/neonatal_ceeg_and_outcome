#!/usr/bin/env Rscript
#
# 07_shap_stability.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Stability of SHAP feature importance.
#   (a) across the 50 outer CV models (SHAP computed on each held-out fold,
#       using the model trained on the remaining folds), and
#   (b) across bootstrap refits (n = CFG$n_boot_shap) of the full-cohort model,
#       with hyperparameters re-tuned inside every bootstrap sample.
# Reported per feature: mean |SHAP|, median rank and 95% rank interval, and
# the proportion of models in which the feature ranks in the top 3.
# ==============================================================================

if (!exists("CFG")) source(file.path("R", "00_config.R"))
if (!exists("run_cv")) source(file.path(CFG$root, "R", "lib_ml.R"))
dat  <- readRDS(file.path(CFG$dir_derived, "analysis_data.rds"))
main <- readRDS(file.path(CFG$dir_derived, "ml_main.rds"))

future::plan(future::multisession, workers = CFG$n_workers)

rank_summary <- function(shap_long, source) {
  shap_long |>
    group_by(outcome, model, unit) |>
    mutate(rank = rank(-mean_abs_shap, ties.method = "average")) |>
    group_by(outcome, model, variable) |>
    summarise(source = source, n_models = n(),
              mean_abs_shap = mean(mean_abs_shap),
              rank_median = median(rank), rank_lo = quantile(rank, .025), rank_hi = quantile(rank, .975),
              top3 = mean(rank <= 3), .groups = "drop")
}

# ---- (a) across outer CV folds -------------------------------------------------
cv_shap <- main$jobs |> filter(algo == "lgb", model != "Clinical") |>
  mutate(shap = map(fit, "shap")) |> select(outcome, model, shap) |> unnest(shap) |>
  mutate(unit = paste(rep, fold))
stab_cv <- rank_summary(cv_shap, "CV folds (50 models)")

# ---- (b) bootstrap refits -----------------------------------------------------
boot_jobs <- expand_grid(outcome = c("neo_death", "epilepsy"), model = c("EEG", "Clinical + EEG"))
feats <- list("EEG" = names(CFG$eeg_all), "Clinical + EEG" = c(names(CFG$clin), names(CFG$eeg_all)))
boot_shap <- pmap_dfr(boot_jobs, function(outcome, model) {
  d <- filter(dat, !is.na(.data[[outcome]]))
  enc <- encode(d, feats[[model]])
  furrr::future_map_dfr(seq_len(CFG$n_boot_shap), function(b) {
    set.seed(CFG$seed + 300000L + b)
    # stratified bootstrap keeps the number of events fixed
    idx <- unlist(lapply(split(seq_len(nrow(d)), d[[outcome]]), function(i) sample(i, replace = TRUE)))
    rec <- prep_learn(enc$X[idx, ])
    Xb <- prep_apply(rec, enc$X[idx, ]); Xall <- prep_apply(rec, enc$X)
    fit <- fit_lgb(Xb, d[[outcome]][idx], seed = CFG$seed + b)
    s <- shap_fit(fit, Xall, enc$map)
    tibble(outcome = outcome, model = model, unit = b, variable = colnames(s),
           mean_abs_shap = colMeans(abs(s)))
  }, .options = furrr::furrr_options(seed = TRUE, packages = c("tidyverse", "lightgbm", "glmnet")))
})
stab_boot <- rank_summary(boot_shap, sprintf("Bootstrap refits (%d models)", CFG$n_boot_shap))

stab <- bind_rows(stab_cv, stab_boot) |>
  mutate(label = coalesce(c(CFG$eeg_all, CFG$clin)[variable], variable)) |>
  arrange(outcome, model, source, rank_median)
write_tsv(stab, file.path(CFG$dir_tab, "tab_shap_stability.tsv"))

future::plan(future::sequential)
saveRDS(list(stab = stab, cv_shap = cv_shap, boot_shap = boot_shap), file.path(CFG$dir_derived, "shap_stability.rds"))
message("[07] SHAP stability done")
