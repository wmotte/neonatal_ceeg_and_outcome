#!/usr/bin/env Rscript
#
# 03_univariable.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# For every outcome x cEEG feature:
#   (a) crude OR, continuous features per 1 SD (standardised), binary as is
#   (b) Benjamini-Hochberg FDR across the 14 features within each outcome
#   (c) OR adjusted for PMA at cEEG
#   (d) OR adjusted for PMA at cEEG and vigilance state
# Firth-penalised logistic regression (logistf) with profile-likelihood CIs,
# because several cells are (quasi-)separated (e.g. recorded seizures vs
# cognition gave an OR of 6.96e4 in the original Table 3).
# Available-case analysis: no imputation is needed for one-feature models.
# For comparison, (e) repeats (a) on the originally imputed features.
# ==============================================================================

if (!exists("CFG")) source(file.path("R", "00_config.R"))
suppressPackageStartupMessages({ library(logistf) })
dat <- readRDS(file.path(CFG$dir_derived, "analysis_data.rds"))

std <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

fit_or <- function(df, outcome, feature, adjust = character()) {
  d <- df |> select(all_of(c(outcome, feature, adjust))) |> drop_na()
  if (!all(d[[feature]] %in% c(0, 1))) d[[feature]] <- std(d[[feature]])
  f <- reformulate(c(feature, adjust), response = outcome)
  m <- logistf(f, data = d, control = logistf.control(maxit = 200),
               plcontrol = logistf.control(maxit = 400))
  i <- which(names(coef(m)) == feature)
  tibble(n = nrow(d), events = sum(d[[outcome]]),
         or = exp(coef(m)[i]), lo = exp(m$ci.lower[i]), hi = exp(m$ci.upper[i]),
         p = m$prob[i])
}

grid <- expand_grid(outcome = names(CFG$outcomes), feature = names(CFG$eeg_all))
uni <- grid |>
  mutate(crude   = map2(outcome, feature, ~ fit_or(dat, .x, .y)),
         adj_pma = map2(outcome, feature, ~ fit_or(dat, .x, .y, "pma")),
         adj_pma_vig = map2(outcome, feature, ~ fit_or(dat, .x, .y, c("pma", "vigilance")))) |>
  pivot_longer(c(crude, adj_pma, adj_pma_vig), names_to = "model", values_to = "res") |>
  unnest(res) |>
  group_by(outcome, model) |>
  mutate(q_fdr = p.adjust(p, method = "BH")) |>
  ungroup() |>
  mutate(label = CFG$eeg_all[feature],
         scale = if_else(feature %in% names(CFG$eeg_visual), "binary", "per SD"))

# (e) crude on the originally imputed features, mirroring the submitted Table 3
imp <- dat |>
  mutate(across(all_of(names(CFG$eeg_quant)), ~ replace_na(.x, median(.x, na.rm = TRUE))),
         bgp_abnormal = replace_na(bgp_abnormal, 0L), asymmetry = replace_na(asymmetry, 0L))
uni_imp <- grid |>
  mutate(res = map2(outcome, feature, ~ fit_or(imp, .x, .y))) |> unnest(res) |>
  group_by(outcome) |> mutate(q_fdr = p.adjust(p, "BH")) |> ungroup() |>
  mutate(model = "crude_imputed", label = CFG$eeg_all[feature],
         scale = if_else(feature %in% names(CFG$eeg_visual), "binary", "per SD"))

uni_all <- bind_rows(uni, uni_imp)
write_tsv(uni_all, file.path(CFG$dir_tab, "tab_univariable_or.tsv"))

# Compact manuscript-style table (crude, FDR, adjusted)
fmt_or <- function(or, lo, hi) sprintf("%.2f (%.2f-%.2f)", or, lo, hi)
fmt_p  <- function(p) if_else(p < 0.001, "<0.001", sprintf("%.3f", p))
tab3 <- uni |>
  mutate(cell = paste0(fmt_or(or, lo, hi), "; p=", fmt_p(p), "; q=", fmt_p(q_fdr))) |>
  select(outcome, feature, label, scale, model, cell) |>
  pivot_wider(names_from = c(outcome, model), values_from = cell)
write_tsv(tab3, file.path(CFG$dir_tab, "tab3_revised_univariable.tsv"))

saveRDS(uni_all, file.path(CFG$dir_derived, "univariable.rds"))
message("[03] univariable ORs written")
