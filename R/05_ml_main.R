#!/usr/bin/env Rscript
#
# 05_ml_main.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Prediction analysis.
# For each outcome, three models replace the "propensity-adjusted" model:
#   Clinical        : GA, PMA at cEEG, neonatal seizures, brain injury, medication
#   EEG             : the 14 cEEG features, unadjusted
#   Clinical + EEG  : both
# Two algorithms: LightGBM (primary) and ridge logistic regression.
# Validation: stratified 5-fold CV x 10 repeats, all preprocessing and tuning
# inside the training folds (lib_ml.R).
# Outputs: tab4_revised_performance.tsv, tab_delta_auc.tsv, tab_permutation.tsv,
#          data/derived/ml_main.rds
# ==============================================================================

if (!exists("CFG")) source(file.path("R", "00_config.R"))
if (!exists("run_cv")) source(file.path(CFG$root, "R", "lib_ml.R"))
dat <- readRDS(file.path(CFG$dir_derived, "analysis_data.rds"))

future::plan(future::multisession, workers = CFG$n_workers)

feature_sets <- list(
  "Clinical"       = names(CFG$clin),
  "EEG"            = names(CFG$eeg_all),
  "Clinical + EEG" = c(names(CFG$clin), names(CFG$eeg_all)))

jobs <- expand_grid(outcome = names(CFG$outcomes), model = names(feature_sets),
                    algo = c("lgb", "ridge"))

# Same seed for every job -> identical outer splits within an outcome, which
# makes the paired AUC comparisons valid.
fits <- furrr::future_pmap(jobs, function(outcome, model, algo) {
  d <- filter(dat, !is.na(.data[[outcome]]))
  run_cv(d, outcome, list(features = feature_sets[[model]], algo = algo))
}, .options = furrr::furrr_options(seed = TRUE, packages = c("tidyverse", "lightgbm", "glmnet")))
jobs$fit <- fits

perf <- jobs |>
  mutate(n = map_int(fit, ~ n_distinct(.x$oof$row)),
         events = map_int(fit, ~ sum(distinct(.x$oof, row, y)$y)),
         summ = map(fit, ~ summarise_cv(.x$oof))) |>
  select(-fit) |> unnest(summ)
write_tsv(perf, file.path(CFG$dir_tab, "tab4_revised_performance.tsv"))

# ---- Incremental value of cEEG: Delta AUC (paired bootstrap) -----------------
get_oof <- function(o, m, a) jobs$fit[[which(jobs$outcome == o & jobs$model == m & jobs$algo == a)]]$oof
dauc <- expand_grid(outcome = names(CFG$outcomes), algo = c("lgb", "ridge")) |>
  mutate(comb_vs_clin = map2(outcome, algo, ~ delta_auc(get_oof(.x, "Clinical + EEG", .y), get_oof(.x, "Clinical", .y))),
         eeg_vs_clin  = map2(outcome, algo, ~ delta_auc(get_oof(.x, "EEG", .y), get_oof(.x, "Clinical", .y)))) |>
  pivot_longer(c(comb_vs_clin, eeg_vs_clin), names_to = "comparison", values_to = "res") |>
  unnest(res)
write_tsv(dauc, file.path(CFG$dir_tab, "tab_delta_auc.tsv"))

# ---- Permutation tests (LightGBM; full pipeline re-run every permutation) -----
message("[05] permutation tests: ", CFG$n_perm, " permutations x 6 models")
perm_jobs <- bind_rows(
  tibble(outcome = names(CFG$outcomes), model = "EEG", type = "label"),
  tibble(outcome = names(CFG$outcomes), model = "Clinical + EEG", type = "eeg"))
perm_jobs$null <- pmap(perm_jobs, function(outcome, model, type) {
  d <- filter(dat, !is.na(.data[[outcome]]))
  perm_test(d, outcome, list(features = feature_sets[[model]], algo = "lgb"), type = type)
})
perm <- perm_jobs |>
  left_join(perf |> filter(algo == "lgb", metric == "auc") |> select(outcome, model, auc = estimate),
            by = c("outcome", "model")) |>
  mutate(null_mean = map_dbl(null, mean), null_q95 = map_dbl(null, ~ quantile(.x, .95)),
         p_perm = map2_dbl(null, auc, ~ (1 + sum(.x >= .y)) / (1 + length(.x))),
         null_hypothesis = if_else(type == "label", "no association (outcome permuted)",
                                   "no added value beyond clinical (cEEG block permuted)"))
write_tsv(select(perm, -null), file.path(CFG$dir_tab, "tab_permutation.tsv"))

future::plan(future::sequential)
saveRDS(list(jobs = jobs, perf = perf, dauc = dauc, perm = perm), file.path(CFG$dir_derived, "ml_main.rds"))
message("[05] main ML analysis done")
