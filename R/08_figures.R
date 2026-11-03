#!/usr/bin/env Rscript
#
# 08_figures.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Revision figures (ggplot2, house style in plot_theme.R).
#   fig_R1_roc                : ROC curves, clinical / EEG / combined (LightGBM)
#   fig_R2_calibration        : calibration of the LightGBM models
#   fig_R3_shap_stability     : SHAP rank stability (CV folds and bootstrap)
#   fig_R4_forest_or          : standardised ORs, crude vs PMA + vigilance adjusted
#   fig_R5_sensitivity        : AUC across sensitivity analyses
#   figS_legacy_asymmetry     : why asymmetry dominated the original SHAP plot
#   figS_qeeg_pma             : quantitative cEEG features versus PMA at cEEG
#   figS_permutation          : permutation null distributions
# ==============================================================================

if (!exists("CFG")) source(file.path("R", "00_config.R"))
if (!exists("theme_neo")) source(file.path(CFG$root, "R", "plot_theme.R"))
if (!exists("auc_fast")) source(file.path(CFG$root, "R", "lib_ml.R"))
suppressPackageStartupMessages({ library(scales) })

dat    <- readRDS(file.path(CFG$dir_derived, "analysis_data.rds"))
main   <- readRDS(file.path(CFG$dir_derived, "ml_main.rds"))
sensr  <- readRDS(file.path(CFG$dir_derived, "ml_sensitivity.rds"))
shp    <- readRDS(file.path(CFG$dir_derived, "shap_stability.rds"))
uni    <- readRDS(file.path(CFG$dir_derived, "univariable.rds"))
legacy <- readRDS(file.path(CFG$dir_derived, "legacy.rds"))

theme_set(theme_neo())
out_lab <- function(x) factor(CFG$outcomes[x], levels = CFG$outcomes)
mod_lev <- names(model_colours)

# subject-level OOF probability (mean over the 10 CV repeats)
oof_avg <- main$jobs |> filter(algo == "lgb") |>
  mutate(o = map(fit, ~ .x$oof |> group_by(row, y) |> summarise(prob = mean(prob), .groups = "drop"))) |>
  select(outcome, model, o) |> unnest(o)

# ---- R1 ROC --------------------------------------------------------------------
roc_df <- oof_avg |> group_by(outcome, model) |>
  group_modify(function(d, k) {
    th <- sort(unique(c(-Inf, d$prob, Inf)), decreasing = TRUE)
    tibble(fpr = map_dbl(th, ~ mean(d$prob[d$y == 0] >= .x)),
           tpr = map_dbl(th, ~ mean(d$prob[d$y == 1] >= .x)))
  }) |> ungroup()
auc_lab <- main$perf |> filter(algo == "lgb", metric == "auc") |>
  mutate(txt = sprintf("%s: %.2f (%.2f-%.2f)", model, estimate, lo, hi)) |>
  group_by(outcome) |> summarise(txt = paste(txt, collapse = "\n"), .groups = "drop")
fig1 <- ggplot(roc_df, aes(fpr, tpr, colour = factor(model, mod_lev))) +
  geom_abline(colour = bc_palette["grey"], linetype = 2) +
  geom_step(linewidth = .8, direction = "vh") +
  geom_text(data = auc_lab, aes(1, 0, label = txt), inherit.aes = FALSE,
            hjust = 1, vjust = 0, size = 2.8, lineheight = .95) +
  facet_wrap(~ out_lab(outcome)) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_x_continuous(breaks = number_ticks(5), labels = percent_format(1)) +
  scale_y_continuous(breaks = number_ticks(5), labels = percent_format(1)) +
  coord_equal() +
  labs(x = "1 - specificity", y = "Sensitivity")
save_fig(fig1, "fig_R1_roc", width = 10, height = 4.2)

# ---- R2 Calibration ------------------------------------------------------------
cal_df <- oof_avg |> filter(model != "Clinical", outcome != "cognition") |>
  group_by(outcome, model) |>
  mutate(bin = ntile(prob, 5)) |>
  group_by(outcome, model, bin) |>
  summarise(pred = mean(prob), obs = mean(y), n = n(),
            lo = binom.test(sum(y), n())$conf.int[1], hi = binom.test(sum(y), n())$conf.int[2],
            .groups = "drop")
fig2 <- ggplot(cal_df, aes(pred, obs, colour = factor(model, mod_lev))) +
  geom_abline(colour = bc_palette["grey"], linetype = 2) +
  geom_smooth(data = oof_avg |> filter(model != "Clinical", outcome != "cognition"),
              aes(prob, y), method = "loess", formula = y ~ x, se = FALSE, linewidth = .6, span = 1) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = .3, position = position_dodge(width = .02)) +
  facet_wrap(~ out_lab(outcome)) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = number_ticks(5)) +
  scale_y_continuous(limits = c(0, 1), breaks = number_ticks(5)) +
  coord_equal() +
  labs(x = "Predicted probability (out-of-fold)", y = "Observed proportion (quintiles, 95% CI)")
save_fig(fig2, "fig_R2_calibration", width = 8, height = 4.4)

# ---- R3 SHAP stability ---------------------------------------------------------
st <- shp$stab |> filter(model == "EEG", outcome != "cognition") |>
  mutate(src = if_else(str_detect(source, "CV"), "Across CV folds", "Across bootstrap refits"))
ord <- st |> filter(src == "Across CV folds") |> group_by(variable) |>
  summarise(r = mean(rank_median)) |> arrange(desc(r)) |> pull(variable)
fig3 <- ggplot(st, aes(rank_median, factor(variable, ord, labels = CFG$eeg_all[ord]), colour = src)) +
  geom_pointrange(aes(xmin = rank_lo, xmax = rank_hi), position = position_dodge(width = .6), size = .3) +
  facet_wrap(~ out_lab(outcome)) +
  scale_x_continuous(breaks = c(1, 5, 10, 14)) +
  scale_colour_manual(values = unname(bc_palette[c("orange", "skyblue")]), name = NULL) +
  labs(x = "SHAP importance rank (median, 95% interval; 1 = most important)", y = NULL)
save_fig(fig3, "fig_R3_shap_stability", width = 9, height = 5.2)

# ---- R4 Forest plot of standardised ORs ----------------------------------------
fo <- uni |> filter(model %in% c("crude", "adj_pma_vig")) |>
  mutate(model = factor(model, c("crude", "adj_pma_vig"),
                        c("Crude", "Adjusted for PMA and vigilance state")),
         sig = if_else(model == "Crude" & q_fdr < 0.05, "FDR q < 0.05", "q >= 0.05"),
         feat = factor(feature, rev(names(CFG$eeg_all)), labels = rev(CFG$eeg_all)))
fig4 <- ggplot(fo, aes(or, feat, colour = model)) +
  geom_vline(xintercept = 1, colour = bc_palette["grey"], linetype = 2) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = .6), size = .25) +
  facet_wrap(~ out_lab(outcome)) +
  scale_x_log10(breaks = c(.1, .25, .5, 1, 2, 4, 10, 25), labels = label_number(drop0trailing = TRUE)) +
  coord_cartesian(xlim = c(.03, 80)) +
  scale_colour_manual(values = unname(bc_palette[c("orange", "blue")]), name = NULL) +
  labs(x = "Odds ratio (Firth; continuous features per 1 SD)", y = NULL)
save_fig(fig4, "fig_R4_forest_or", width = 10, height = 5.2)

# ---- R5 Sensitivity analyses ---------------------------------------------------
se <- sensr$sens |> filter(metric == "auc", outcome != "cognition") |>
  mutate(lab = sprintf("%s: %s (n=%d, events=%d)", scenario, label, n, events))
fig5 <- ggplot(se, aes(estimate, fct_rev(factor(scenario)), colour = factor(model, mod_lev))) +
  geom_vline(xintercept = .5, colour = bc_palette["grey"], linetype = 2) +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = .6), size = .25) +
  facet_wrap(~ out_lab(outcome)) +
  scale_colour_manual(values = model_colours, name = NULL) +
  scale_x_continuous(breaks = number_ticks(5)) +
  labs(x = "AUC (95% CI), LightGBM, repeated nested CV", y = NULL,
       caption = paste(unique(paste0(se$scenario, " = ", se$label)), collapse = "; ")) +
  theme(plot.caption = element_text(hjust = 0, size = 7))
save_fig(fig5, "fig_R5_sensitivity", width = 10, height = 5)

# ---- S: legacy asymmetry artefact -------------------------------------------------
sy <- legacy$sym |> filter(outcome != "cognition") |>
  mutate(Asymmetry = factor(asymmetry, 0:1, c("Symmetric", "Asymmetric")),
         Outcome = factor(y, 0:1, c("No", "Yes")))
figS1 <- ggplot(sy, aes(ps, asym_resid, colour = Outcome, shape = Asymmetry)) +
  geom_point(size = 1.8, alpha = .8) +
  facet_wrap(~ out_lab(outcome)) +
  scale_colour_manual(values = unname(bc_palette[c("skyblue", "vermillion")])) +
  scale_x_continuous(breaks = number_ticks(5)) +
  labs(x = "Prognostic ('propensity') score from the clinical confounders",
       y = "Asymmetry feature after OLS residualisation (z)",
       colour = "Outcome", shape = NULL) +
  theme(legend.title = element_text(size = rel(.9)))
save_fig(figS1, "figS_legacy_asymmetry", width = 8.5, height = 4.2)

# ---- S: qEEG vs PMA --------------------------------------------------------------
qp <- dat |> select(record, pma, vigilance, all_of(names(CFG$eeg_quant))) |>
  pivot_longer(all_of(names(CFG$eeg_quant)), names_to = "feature", values_to = "value") |>
  filter(!is.na(value), !is.na(vigilance)) |>
  mutate(feature = factor(feature, names(CFG$eeg_quant), CFG$eeg_quant))
figS2 <- ggplot(qp, aes(pma, value)) +
  geom_point(aes(colour = vigilance), size = 1.2, alpha = .8) +
  geom_smooth(method = "loess", formula = y ~ x, colour = "black", linewidth = .6, se = FALSE) +
  facet_wrap(~ feature, scales = "free_y", ncol = 5) +
  scale_colour_manual(values = unname(bc_palette[c("orange", "skyblue", "green", "vermillion")]), name = NULL) +
  scale_x_continuous(breaks = number_ticks(4)) +
  labs(x = "PMA at cEEG (weeks)", y = NULL)
save_fig(figS2, "figS_qeeg_pma", width = 11, height = 5)

# ---- S: permutation nulls ---------------------------------------------------------
pn <- main$perm |> select(outcome, model, type, auc, p_perm, null) |> unnest(null) |>
  mutate(panel = sprintf("%s\n%s", CFG$outcomes[outcome],
                         if_else(type == "label", "EEG model, outcome permuted", "Clinical + EEG, cEEG permuted")))
pl <- distinct(pn, panel, auc, p_perm)
figS3 <- ggplot(pn, aes(null)) +
  geom_histogram(bins = 40, fill = bc_palette["grey"], colour = "white") +
  geom_vline(data = pl, aes(xintercept = auc), colour = bc_palette["vermillion"], linewidth = .8) +
  geom_text(data = pl, aes(x = auc, y = Inf, label = sprintf("p = %.3f", p_perm)),
            hjust = 1.1, vjust = 1.5, size = 3) +
  facet_wrap(~ panel, ncol = 3) +
  labs(x = "Cross-validated AUC under the null (full pipeline per permutation)", y = "Count")
save_fig(figS3, "figS_permutation", width = 10, height = 5.5)

message("[08] figures written to ", CFG$dir_fig)
