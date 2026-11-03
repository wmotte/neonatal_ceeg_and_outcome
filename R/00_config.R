#!/usr/bin/env Rscript
#
# 00_config.R
# Central configuration.
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Every path, seed, variable set, and analysis setting lives here, so that the
# numbered scripts (01-08) contain no hard-coded magic numbers.
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ---- Paths -------------------------------------------------------------------
# Resolve the project root by searching upward for R/00_config.R, so the
# pipeline runs from any working directory (source() or Rscript).
CFG <- list()
.find_root <- function() {
  env <- Sys.getenv("NEOEEG_ROOT", "")
  if (nzchar(env) && file.exists(file.path(env, "R", "00_config.R"))) return(normalizePath(env))
  args <- commandArgs(FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  start <- if (length(fa)) dirname(sub("^--file=", "", fa[1])) else getwd()
  d <- normalizePath(start, mustWork = FALSE)
  for (i in 1:6) {
    if (file.exists(file.path(d, "R", "00_config.R"))) return(d)
    parent <- dirname(d)
    if (parent == d) break
    d <- parent
  }
  normalizePath(getwd(), mustWork = FALSE)
}
CFG$root        <- .find_root()
CFG$dir_raw     <- file.path(CFG$root, "data", "raw")
CFG$dir_legacy  <- file.path(CFG$root, "data", "legacy")
CFG$dir_derived <- file.path(CFG$root, "data", "derived")
CFG$dir_tab     <- file.path(CFG$root, "output", "tables")
CFG$dir_fig     <- file.path(CFG$root, "output", "figures")
CFG$dir_diag    <- file.path(CFG$root, "output", "diagnostics")
CFG$sav         <- file.path(CFG$dir_raw, "20260330_Database_neonatal_EEG.sav")
CFG$legacy_csv  <- file.path(CFG$dir_legacy, "Dataset_20231124.csv")

for (d in c(CFG$dir_derived, CFG$dir_tab, CFG$dir_fig, CFG$dir_diag)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- Reproducibility ---------------------------------------------------------
CFG$seed      <- 20261103L
CFG$n_workers <- max(1L, min(8L, parallel::detectCores() - 2L))

# ---- Outcomes ----------------------------------------------------------------
CFG$outcomes <- c(neo_death = "Neonatal death",
                  epilepsy  = "Epilepsy",
                  cognition = "Cognitive impairment")

# ---- Feature sets ------------------------------------------------------------
# The 14 cEEG features of the original analysis (4 visual + 10 quantitative).
CFG$eeg_visual <- c(bgp_abnormal = "Abnormal background pattern",
                    asymmetry    = "Asymmetry",
                    ieds         = "Interictal epileptiform discharges",
                    rec_seiz     = "Recorded seizures")
CFG$eeg_quant  <- c(power_delta = "Spectral power, delta",
                    power_theta = "Spectral power, theta",
                    corr_delta  = "Correlation, delta",
                    corr_theta  = "Correlation, theta",
                    coh_delta   = "Coherence, delta",
                    coh_theta   = "Coherence, theta",
                    reeg_delta  = "Range EEG median, delta",
                    reeg_theta  = "Range EEG median, theta",
                    mse_auc     = "Multiscale entropy, AUC",
                    mse_max     = "Multiscale entropy, maximum")
CFG$eeg_all    <- c(CFG$eeg_visual, CFG$eeg_quant)

# Clinical predictors available at the time of the cEEG. The four original
# "confounders" plus PMA at cEEG (reviewer 3, comment 3). GA at birth and PMA
# together also encode postnatal age at recording.
CFG$clin       <- c(ga          = "Gestational age at birth",
                    pma         = "PMA at cEEG",
                    seizure_cat = "Neonatal seizures (aetiology)",
                    brain_injury = "Moderate-severe brain injury",
                    medication  = "Medication affecting EEG")
# Vigilance state is examined as a sensitivity covariate only (see 06_*).
CFG$vigilance  <- c(vigilance = "Vigilance state")

# ---- Cross-validation and tuning ---------------------------------------------
CFG$cv_folds    <- 5L     # outer folds (stratified)
CFG$cv_repeats  <- 10L    # outer repeats
CFG$inner_folds <- 5L     # inner folds for hyperparameter tuning
# LightGBM grid; deliberately small (n = 65-127). Tuned on inner-fold log-loss.
CFG$lgb_grid <- tidyr::expand_grid(max_depth = c(1L, 2L),
                                   nrounds = c(50L, 100L, 200L),
                                   min_data_in_leaf = c(5L, 10L)) |>
  dplyr::mutate(learning_rate = 0.05, num_leaves = 2L^max_depth)
CFG$n_boot      <- 2000L  # bootstrap resamples for metric CIs
CFG$n_perm      <- 1000L  # permutations (full pipeline re-run each time)
CFG$n_boot_shap <- 200L   # bootstrap refits for SHAP stability

# ---- Sensitivity thresholds --------------------------------------------------
CFG$min_channels <- 12L   # exclude recordings with < 12 of 18 channels retained
CFG$min_pma_term <- 37    # PMA >= 37 weeks at cEEG subgroup
CFG$min_fu_epi_y <- 5     # epilepsy-free children need >= 5 years follow-up

# ---- Deaths: cEEG documented as part of a redirection-of-care decision -------
# Coded by WO from the free-text field `cause_death` (SPSS truncates it at 255
# characters). To be verified by LN against the patient records.
CFG$death_eeg_cited <- c("REC00683", "REC00725", "REC00829", "REC00889",
                         "REC01001", "REC01087", "REC01099", "REC01192",
                         "REC01291", "REC01335", "REC01349", "REC01609",
                         "REC02391", "REC02398")
CFG$death_eeg_uncertain <- c("REC01255", "REC02213", "REC02214", "REC02403")
