#!/usr/bin/env Rscript
#
# 06_ml_sensitivity.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Sensitivity analyses (LightGBM; EEG and Clinical + EEG models):
#   S1  + vigilance state as clinical covariate                  
#   S2  quantitative cEEG features normalised for PMA within
#       the training folds                                       
#   S3  PMA >= 37 weeks at cEEG                                  
#   S4  excluding altered state of consciousness                
#   S5  complete cases for all quantitative cEEG features         
#   S6  excluding recordings with < 12 retained channels         
#   S7  neonatal death: excluding deaths where the cEEG was cited
#       (or possibly cited) in the redirection-of-care decision   
#   S8  epilepsy: epilepsy-free children need >= 5 y follow-up    
# Plus discrimination of the main models separately for deaths after
# redirection of care and for other deaths, using the main out-of-fold
# predictions, and correlation of predictions with channel count.
# ==============================================================================

if (!exists("CFG")) source(file.path("R", "00_config.R"))
if (!exists("run_cv")) source(file.path(CFG$root, "R", "lib_ml.R"))
dat  <- readRDS(file.path(CFG$dir_derived, "analysis_data.rds"))
main <- readRDS(file.path(CFG$dir_derived, "ml_main.rds"))

future::plan(future::multisession, workers = CFG$n_workers)

clin <- names(CFG$clin); eeg <- names(CFG$eeg_all)
scen <- list(
  S1 = list(label = "+ vigilance state", filter = function(d, o) d,
            feat = list("EEG" = c(eeg, "vigilance"), "Clinical + EEG" = c(clin, "vigilance", eeg))),
  S2 = list(label = "qEEG normalised for PMA (within folds)", filter = function(d, o) d, pma_norm = TRUE,
            feat = list("EEG" = c(eeg), "Clinical + EEG" = c(clin, eeg))),
  S3 = list(label = sprintf("PMA >= %d weeks at cEEG", CFG$min_pma_term),
            filter = function(d, o) filter(d, pma >= CFG$min_pma_term)),
  S4 = list(label = "Excluding altered consciousness", filter = function(d, o) filter(d, altered_consc == 0)),
  S5 = list(label = "Complete quantitative cEEG", filter = function(d, o) filter(d, if_all(all_of(names(CFG$eeg_quant)), ~ !is.na(.x)))),
  S6 = list(label = sprintf("Excluding < %d channels", CFG$min_channels),
            filter = function(d, o) filter(d, is.na(n_channels) | n_channels >= CFG$min_channels)),
  S7 = list(label = "Excluding deaths with cEEG (possibly) cited in redirection decision", only = "neo_death",
            filter = function(d, o) filter(d, !(death_eeg_cited %in% c("cEEG cited", "Uncertain")))),
  S8 = list(label = sprintf("Epilepsy-free children with >= %d y follow-up", CFG$min_fu_epi_y), only = "epilepsy",
            filter = function(d, o) filter(d, epilepsy == 1 | fu_years >= CFG$min_fu_epi_y))
)
default_feat <- list("EEG" = eeg, "Clinical + EEG" = c(clin, eeg))

jobs <- imap_dfr(scen, function(s, id) {
  outs <- s$only %||% names(CFG$outcomes)
  ff <- s$feat %||% default_feat
  expand_grid(scenario = id, outcome = outs, model = names(ff)) |>
    mutate(label = s$label)
})

fits <- furrr::future_pmap(select(jobs, scenario, outcome, model), function(scenario, outcome, model) {
  s <- scen[[scenario]]
  ff <- s$feat %||% default_feat
  d <- s$filter(filter(dat, !is.na(.data[[outcome]])), outcome)
  d <- droplevels(d)
  # drop constant columns (e.g. vigilance levels absent in a subgroup)
  r <- run_cv(d, outcome, list(features = ff[[model]], algo = "lgb", pma_normalise = isTRUE(s$pma_norm)),
              keep_shap = FALSE)
  list(oof = r$oof, n = nrow(d), events = sum(d[[outcome]]))
}, .options = furrr::furrr_options(seed = TRUE, packages = c("tidyverse", "lightgbm", "glmnet")))
jobs$n <- map_int(fits, "n"); jobs$events <- map_int(fits, "events")
jobs$summ <- map(fits, ~ summarise_cv(.x$oof, n_boot = 1000L))
sens <- unnest(jobs, summ) |> filter(metric %in% c("auc", "bacc", "sens", "spec", "brier"))

# reference rows from the main analysis
ref <- main$perf |> filter(algo == "lgb", model %in% c("EEG", "Clinical + EEG"),
                           metric %in% c("auc", "bacc", "sens", "spec", "brier")) |>
  transmute(scenario = "S0", outcome, model, label = "Main analysis", n, events, metric, estimate, lo, hi)
sens <- bind_rows(ref, sens)
write_tsv(sens, file.path(CFG$dir_tab, "tab_sensitivity.tsv"))

# ---- R3.2 Redirection of care: discrimination by type of death ----------------
oof_main <- function(m) {
  main$jobs |> filter(outcome == "neo_death", model == m, algo == "lgb") |> pull(fit) |> pluck(1, "oof") |>
    group_by(row, record, y) |> summarise(prob = mean(prob), thr = mean(thr), .groups = "drop") |>
    left_join(select(dat, record, death_redirection, death_eeg_cited), by = "record")
}
redir <- map_dfr(c("EEG", "Clinical + EEG"), function(m) {
  o <- oof_main(m)
  surv <- filter(o, y == 0)
  grp <- list("Deaths after redirection of care" = filter(o, y == 1, death_redirection == 1),
              "Other neonatal deaths" = filter(o, y == 1, death_redirection == 0),
              "Deaths, cEEG cited in decision" = filter(o, y == 1, death_eeg_cited == "cEEG cited"),
              "Deaths, cEEG not (clearly) cited" = filter(o, y == 1, death_eeg_cited != "cEEG cited"))
  imap_dfr(grp, function(g, nm) {
    yy <- c(rep(1, nrow(g)), rep(0, nrow(surv))); pp <- c(g$prob, surv$prob)
    set.seed(CFG$seed)
    bt <- replicate(CFG$n_boot, {
      a <- sample(nrow(g), replace = TRUE); b <- sample(nrow(surv), replace = TRUE)
      auc_fast(c(rep(1, nrow(g)), rep(0, nrow(surv))), c(g$prob[a], surv$prob[b]))
    })
    tibble(model = m, group = nm, n_deaths = nrow(g), n_survivors = nrow(surv),
           auc_vs_survivors = auc_fast(yy, pp), lo = quantile(bt, .025), hi = quantile(bt, .975),
           sensitivity = mean(g$prob >= g$thr), median_prob = median(g$prob))
  })
})
write_tsv(redir, file.path(CFG$dir_tab, "tab_death_by_redirection.tsv"))

# ---- R3.7 Model output vs number of retained channels ------------------------
chan <- main$jobs |> filter(algo == "lgb", model %in% c("EEG", "Clinical + EEG")) |>
  mutate(res = pmap(list(fit, outcome), function(f, o) {
    p <- f$oof |> group_by(record, y) |> summarise(prob = mean(prob), .groups = "drop") |>
      left_join(select(dat, record, n_channels), by = "record") |> filter(!is.na(n_channels))
    # partial association: residualise prob on outcome first
    pr <- resid(lm(logit(clip01(prob)) ~ y, data = p))
    ct <- suppressWarnings(cor.test(pr, p$n_channels, method = "spearman"))
    tibble(n = nrow(p), rho_channels_given_outcome = unname(ct$estimate), p = ct$p.value)
  })) |> select(outcome, model, res) |> unnest(res)
write_tsv(chan, file.path(CFG$dir_tab, "tab_prediction_vs_channels.tsv"))

future::plan(future::sequential)
saveRDS(list(sens = sens, redir = redir, chan = chan), file.path(CFG$dir_derived, "ml_sensitivity.rds"))
message("[06] sensitivity analyses done")
