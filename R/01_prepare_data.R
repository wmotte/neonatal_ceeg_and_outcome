#!/usr/bin/env Rscript
#
# 01_prepare_data.R
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
#
# Reads the SPSS database (20260330), recodes the variables used in the
# revision, and writes data/derived/analysis_data.rds.
#
# Important: we start from the RAW (non-imputed) cEEG features. The *_new and
# *_interpolation columns in the SPSS file were imputed on the full cohort; all
# imputation is now done inside the training folds (see lib_ml.R).
# ==============================================================================

if (!exists("CFG")) source(file.path("R", "00_config.R"))
suppressPackageStartupMessages({ library(haven) })

raw <- read_sav(CFG$sav)

na_codes <- function(x) { x <- zap_labels(x); x[x %in% c(555, 777, 999)] <- NA; x }

dat <- raw |>
  transmute(
    record,
    sex          = factor(zap_labels(sex), 1:2, c("male", "female")),
    ga           = ga,                          # weeks
    pma          = pma,                         # weeks, at the included cEEG
    pna_days     = (pma - ga) * 7,              # postnatal age at cEEG
    ga_cat       = factor(zap_labels(ga_cate), 0:2,
                          c("Term", "Moderate-late preterm", "Extremely-very preterm")),
    seizure_cat  = factor(zap_labels(neo_seizure_new), 0:3,
                          c("None", "Hypoxic/structural", "Infectious", "Metabolic/genetic/unknown")),
    brain_injury = as.integer(zap_labels(struct_brain_abn_neo_new2) == 2),
    medication   = as.integer(zap_labels(medication_EEG) == 1),
    intubation   = as.integer(zap_labels(intubation) == 1),
    vigilance    = factor(na_codes(vigilance_state), 1:4,
                          c("Wake", "Sleep", "Wake/sleep", "Altered consciousness")),
    altered_consc = as.integer(vigilance == "Altered consciousness"),
    # outcomes
    neo_death    = as.integer(zap_labels(neo_death)),
    epilepsy     = as.integer(zap_labels(epilepsy_new)),
    cognition    = as.integer(zap_labels(cognition_new)),
    death_any    = as.integer(zap_labels(death)),
    # visual cEEG features (raw, so 'not assessable' -> NA)
    bgp_abnormal = as.integer(na_codes(bgp) == 0),
    asymmetry    = as.integer(na_codes(symmetry) == 0),
    ieds         = as.integer(zap_labels(ieds)),
    rec_seiz     = as.integer(zap_labels(rec_seiz)),
    # quantitative cEEG features (raw, NA where not computable)
    power_delta = spectral_power_1, power_theta = spectral_power_2,
    corr_delta  = connectivity_corr_1, corr_theta = connectivity_corr_2,
    coh_delta   = connectivity_coh_mean_1, coh_theta = connectivity_coh_mean_2,
    reeg_delta  = rEEG_median_1, reeg_theta = rEEG_median_2,
    mse_auc     = MSE_auc, mse_max = MSE_max,
    # recording quality
    n_channels  = resid_ch_num,
    eeg_minutes = resid_avg_minutes,
    # follow-up and death timing
    date_eeg    = date_included_neo_eeg,
    date_death  = date_death,
    eeg_to_death_days = as.numeric(date_death - date_included_neo_eeg),
    fu_years    = number_days_wo_neo_death / 365.25,   # birth -> last FU (neonatal survivors)
    followup_status = as_factor(followup_status),
    cause_death,
    # legacy prognostic scores computed in SPSS on the full cohort
    ps_neo_death = PS_neo_death, ps_epilepsy = PS_epilepsy_new, ps_cognition = PS_cognition_new
  ) |>
  mutate(
    death_redirection = if_else(neo_death == 1L,
                                as.integer(str_detect(cause_death, "^redirection of care")), NA_integer_),
    death_eeg_cited   = case_when(neo_death != 1L ~ NA_character_,
                                  record %in% CFG$death_eeg_cited ~ "cEEG cited",
                                  record %in% CFG$death_eeg_uncertain ~ "Uncertain",
                                  TRUE ~ "Not cited"),
    qeeg_missing = as.integer(is.na(power_delta)),
    conn_missing = as.integer(is.na(corr_delta))
  ) |>
  select(-cause_death)

# ---- Sanity checks against the manuscript ------------------------------------
stopifnot(nrow(dat) == 127,
          sum(dat$neo_death) == 33,
          sum(!is.na(dat$epilepsy)) == 68, sum(dat$epilepsy, na.rm = TRUE) == 28,
          sum(!is.na(dat$cognition)) == 65, sum(dat$cognition, na.rm = TRUE) == 36,
          sum(dat$death_redirection, na.rm = TRUE) == 25,
          !anyNA(dat$ga_cat), !anyNA(dat$sex), !anyNA(dat$vigilance))

# Optional: age at epilepsy onset (needed for a time-to-event analysis).
# Not present in the 20260330 database; 05_descriptives.R checks again.
CFG$has_epi_onset <- "age_epilepsy_onset_days" %in% names(raw)

saveRDS(dat, file.path(CFG$dir_derived, "analysis_data.rds"))
message("[01] analysis data: ", nrow(dat), " neonates, ", ncol(dat), " variables")
