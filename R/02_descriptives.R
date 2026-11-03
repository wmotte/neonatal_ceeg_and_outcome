#!/usr/bin/env Rscript
#
# 02_descriptives.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Descriptive analyses:
#   R1.1  missingness per variable                      -> tab_missingness.tsv
#   R1.2  neonates with vs without quantitative cEEG     -> tab_qeeg_missing_vs_complete.tsv
#   R3.2  timing of cEEG relative to death; redirection  -> tab_death_timing.tsv
#   R3.3  qEEG vs PMA and vigilance state; overlap of
#         PMA/GA and vigilance/medication                -> tab_qeeg_by_pma_vigilance.tsv,
#                                                           tab_vigilance_medication.tsv
#   R3.4  follow-up by epilepsy status; survivors with
#         vs without epilepsy outcome                    -> tab_epilepsy_followup.tsv,
#                                                           tab_epilepsy_attrition.tsv
#   R3.7  channel availability vs qEEG features          -> tab_qeeg_by_channels.tsv
# ==============================================================================

if (!exists("CFG")) source(file.path("R", "00_config.R"))
dat <- readRDS(file.path(CFG$dir_derived, "analysis_data.rds"))

fmt_med <- function(x, d = 1) {
  x <- x[!is.na(x)]
  sprintf(paste0("%.", d, "f (%.", d, "f-%.", d, "f) [%.", d, "f-%.", d, "f]"),
          median(x), quantile(x, .25), quantile(x, .75), min(x), max(x))
}
fmt_n  <- function(x) sprintf("%d/%d (%.0f%%)", sum(x, na.rm = TRUE), sum(!is.na(x)),
                              100 * mean(x, na.rm = TRUE))

# Compare two groups on a fixed set of characteristics.
compare_groups <- function(df, group, vars_cont, vars_bin, vars_fac) {
  g <- df[[group]]
  cont <- map_dfr(vars_cont, function(v) tibble(
    variable = v, level = "median (IQR) [range]",
    group0 = fmt_med(df[[v]][g == 0]), group1 = fmt_med(df[[v]][g == 1]),
    p = suppressWarnings(wilcox.test(df[[v]] ~ g)$p.value), test = "Mann-Whitney"))
  bin <- map_dfr(vars_bin, function(v) tibble(
    variable = v, level = "yes",
    group0 = fmt_n(df[[v]][g == 0]), group1 = fmt_n(df[[v]][g == 1]),
    p = fisher.test(table(df[[v]], g))$p.value, test = "Fisher"))
  fac <- map_dfr(vars_fac, function(v) {
    tb <- table(df[[v]], g)
    tibble(variable = v, level = rownames(tb),
           group0 = sprintf("%d (%.0f%%)", tb[, 1], 100 * prop.table(tb[, 1])),
           group1 = sprintf("%d (%.0f%%)", tb[, 2], 100 * prop.table(tb[, 2])),
           p = c(fisher.test(tb, simulate.p.value = TRUE, B = 1e5)$p.value, rep(NA, nrow(tb) - 1)),
           test = c("Fisher (simulated)", rep("", nrow(tb) - 1)))
  })
  bind_rows(cont, bin, fac)
}

# ---- R1.1 Missingness per variable -------------------------------------------
miss_vars <- c(names(CFG$clin), "vigilance", names(CFG$eeg_all), "n_channels",
               "neo_death", "epilepsy", "cognition")
tab_miss <- tibble(variable = miss_vars) |>
  mutate(label = coalesce(c(CFG$clin, CFG$vigilance, CFG$eeg_all)[variable], variable),
         n_total = nrow(dat),
         n_missing = map_int(variable, ~ sum(is.na(dat[[.x]]))),
         pct_missing = round(100 * n_missing / n_total, 1))
# outcome missingness among neonatal survivors is the relevant denominator
surv <- filter(dat, neo_death == 0)
tab_miss <- tab_miss |>
  mutate(n_total   = if_else(variable %in% c("epilepsy", "cognition"), nrow(surv), n_total),
         n_missing = if_else(variable %in% c("epilepsy", "cognition"),
                             map_int(variable, ~ sum(is.na(surv[[.x]]))), n_missing),
         pct_missing = round(100 * n_missing / n_total, 1))
write_tsv(tab_miss, file.path(CFG$dir_tab, "tab_missingness.tsv"))

# ---- R1.2 Neonates with vs without quantitative cEEG --------------------------
tab_qmiss <- compare_groups(dat, "qeeg_missing",
  vars_cont = c("ga", "pma", "pna_days"),
  vars_bin  = c("neo_death", "brain_injury", "medication", "intubation", "bgp_abnormal",
                "rec_seiz", "ieds", "epilepsy", "cognition"),
  vars_fac  = c("seizure_cat", "vigilance")) |>
  rename(qeeg_available = group0, qeeg_missing = group1)
write_tsv(tab_qmiss, file.path(CFG$dir_tab, "tab_qeeg_missing_vs_complete.tsv"))

# ---- R3.2 Timing of cEEG relative to death -----------------------------------
deaths <- filter(dat, neo_death == 1) |>
  mutate(group = if_else(death_redirection == 1, "Redirection of care", "Other death"))
death_summ <- function(df) summarise(df, n = n(),
  eeg_to_death = fmt_med(eeg_to_death_days, 0),
  within_7d    = fmt_n(eeg_to_death_days <= 7),
  pma_at_eeg   = fmt_med(pma, 1),
  eeg_cited    = sum(death_eeg_cited == "cEEG cited"),
  eeg_uncertain = sum(death_eeg_cited == "Uncertain"), .groups = "drop")
tab_death <- bind_rows(
  deaths |> death_summ() |> mutate(group = "All neonatal deaths", .before = 1),
  deaths |> group_by(group) |> death_summ())
write_tsv(tab_death, file.path(CFG$dir_tab, "tab_death_timing.tsv"))

# ---- R3.3 qEEG versus PMA, vigilance state and medication ---------------------
q <- names(CFG$eeg_quant)
tab_pma <- map_dfr(q, function(v) {
  x <- dat[[v]]
  ct <- suppressWarnings(cor.test(x, dat$pma, method = "spearman"))
  cg <- suppressWarnings(cor.test(x, dat$ga, method = "spearman"))
  kw <- kruskal.test(x ~ dat$vigilance)
  mw <- suppressWarnings(wilcox.test(x ~ dat$medication))
  tibble(variable = v, label = CFG$eeg_quant[[v]], n = sum(!is.na(x)),
         rho_pma = unname(ct$estimate), p_pma = ct$p.value,
         rho_ga = unname(cg$estimate), p_ga = cg$p.value,
         p_vigilance_kw = kw$p.value, p_medication_mw = mw$p.value)
}) |> mutate(across(starts_with("p_"), ~ p.adjust(.x, "BH"), .names = "q_{.col}"))
write_tsv(tab_pma, file.path(CFG$dir_tab, "tab_qeeg_by_pma_vigilance.tsv"))

tab_vm <- dat |> count(vigilance, medication) |>
  pivot_wider(names_from = medication, values_from = n, names_prefix = "medication_", values_fill = 0)
write_tsv(tab_vm, file.path(CFG$dir_tab, "tab_vigilance_medication.tsv"))
overlap <- tibble(
  stat = c("spearman_ga_pma", "pearson_ga_pma", "median_pna_days", "cramers_v_vigilance_medication",
           "altered_consc_with_medication"),
  value = c(cor(dat$ga, dat$pma, method = "spearman"), cor(dat$ga, dat$pma),
            median(dat$pna_days),
            sqrt(suppressWarnings(chisq.test(table(dat$vigilance, dat$medication)))$statistic / nrow(dat)),
            sum(dat$altered_consc == 1 & dat$medication == 1) / sum(dat$altered_consc == 1)))
write_tsv(overlap, file.path(CFG$dir_tab, "tab_confounder_overlap.tsv"))

# ---- R3.4 Epilepsy: follow-up and attrition -----------------------------------
epi <- filter(dat, !is.na(epilepsy))
tab_epi_fu <- epi |> group_by(epilepsy) |>
  summarise(n = n(), fu_lt_5y = fmt_n(fu_years < CFG$min_fu_epi_y),
            fu_lt_2y = fmt_n(fu_years < 2), fu_years = fmt_med(fu_years, 1), .groups = "drop") |>
  mutate(p_wilcoxon = wilcox.test(fu_years ~ epilepsy, data = epi)$p.value)
write_tsv(tab_epi_fu, file.path(CFG$dir_tab, "tab_epilepsy_followup.tsv"))

surv <- surv |> mutate(epi_missing = as.integer(is.na(epilepsy)))
tab_epi_attr <- compare_groups(surv, "epi_missing",
  vars_cont = c("ga", "pma", "fu_years"),
  vars_bin  = c("brain_injury", "medication", "bgp_abnormal", "asymmetry", "ieds", "rec_seiz",
                "qeeg_missing"),
  vars_fac  = c("sex", "ga_cat", "seizure_cat", "vigilance")) |>
  rename(epilepsy_outcome_available = group0, epilepsy_outcome_missing = group1)
# quantitative cEEG features too
tab_epi_attr_q <- compare_groups(surv, "epi_missing", vars_cont = names(CFG$eeg_quant),
                                 vars_bin = character(), vars_fac = character()) |>
  rename(epilepsy_outcome_available = group0, epilepsy_outcome_missing = group1)
write_tsv(bind_rows(tab_epi_attr, tab_epi_attr_q), file.path(CFG$dir_tab, "tab_epilepsy_attrition.tsv"))

# Time-to-event analysis needs age at epilepsy onset (not in the database yet).
if (!isTRUE(CFG$has_epi_onset)) {
  message("[02] age at epilepsy onset not available: time-to-event analysis skipped")
}

# ---- R3.7 Channel availability -----------------------------------------------
tab_ch <- map_dfr(q, function(v) {
  ct <- suppressWarnings(cor.test(dat[[v]], dat$n_channels, method = "spearman"))
  tibble(variable = v, label = CFG$eeg_quant[[v]], rho_channels = unname(ct$estimate), p = ct$p.value)
}) |> mutate(q_bh = p.adjust(p, "BH"))
ch_dist <- dat |> filter(!is.na(n_channels)) |> count(n_channels)
write_tsv(tab_ch, file.path(CFG$dir_tab, "tab_qeeg_by_channels.tsv"))
write_tsv(ch_dist, file.path(CFG$dir_tab, "tab_channel_distribution.tsv"))
ch_death <- dat |> filter(!is.na(n_channels)) |>
  summarise(rho_channels_death = cor(n_channels, neo_death, method = "spearman"),
            p = suppressWarnings(wilcox.test(n_channels ~ neo_death)$p.value),
            n_lt_min = sum(n_channels < CFG$min_channels))
write_tsv(ch_death, file.path(CFG$dir_tab, "tab_channels_vs_death.tsv"))

message("[02] descriptives written to ", CFG$dir_tab)
