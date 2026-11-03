#!/usr/bin/env Rscript
#
# 09_report.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Collects all output tables.
# ==============================================================================

if (!exists("CFG")) source(file.path("R", "00_config.R"))

tabs <- list.files(CFG$dir_tab, pattern = "\\.tsv$", full.names = TRUE)
sheets <- set_names(map(tabs, read_tsv, show_col_types = FALSE),
                    str_trunc(str_remove(basename(tabs), "\\.tsv$"), 31, ellipsis = ""))
writexl::write_xlsx(sheets, file.path(CFG$root, "output", "revision_tables.xlsx"))

main  <- readRDS(file.path(CFG$dir_derived, "ml_main.rds"))
sensr <- readRDS(file.path(CFG$dir_derived, "ml_sensitivity.rds"))
f3 <- function(e, l, h) sprintf("%.2f (95%% CI %.2f-%.2f)", e, l, h)

perf_md <- main$perf |>
  filter(metric %in% c("auc", "bacc", "sens", "spec", "ppv", "brier", "cal_intercept", "cal_slope")) |>
  mutate(val = if_else(is.na(lo), sprintf("%.2f", estimate), f3(estimate, lo, hi))) |>
  select(outcome, model, algo, metric, val) |>
  pivot_wider(names_from = metric, values_from = val)
dauc_md <- main$dauc |> mutate(val = sprintf("%+.3f (%.3f to %.3f), p=%.3f", delta_auc, lo, hi, p_boot)) |>
  select(outcome, algo, comparison, val)
perm_md <- main$perm |> select(outcome, model, type, auc, null_mean, null_q95, p_perm)

md_table <- function(df) {
  df <- mutate(df, across(everything(), ~ if (is.numeric(.x)) sprintf("%.3f", .x) else as.character(.x)))
  c(paste("|", paste(names(df), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|"),
    apply(df, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}

lines <- c("# Revision analyses: key numbers", "",
           sprintf("Generated %s by R/09_report.R.", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
           "## Performance (repeated nested CV, mean over 10 repeats; bootstrap 95% CI)", "",
           md_table(perf_md), "",
           "## Incremental value: Delta AUC (paired bootstrap)", "", md_table(dauc_md), "",
           "## Permutation tests (LightGBM)", "", md_table(perm_md), "",
           "## Neonatal death by type of death (main out-of-fold predictions)", "",
           md_table(sensr$redir), "",
           "## Sensitivity analyses (AUC)", "",
           md_table(sensr$sens |> filter(metric == "auc") |>
                      transmute(scenario, label, outcome, model, n, events, auc = f3(estimate, lo, hi))))
writeLines(lines, file.path(CFG$root, "output", "results_summary.md"))
message("[09] report written")
