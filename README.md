# Neonatal cEEG and outcome

## Quick start

```bash
Rscript R/install_requirements.R   # once
Rscript R/run_all.R                # about 20-30 min on 8 cores
```

## Scripts

| Script | Comment |
|---|---|
| `00_config.R` | paths, seed, feature sets, CV settings, thresholds |
| `plot_theme.R` | figure house style (Okabe-Ito, `theme_classic`, `save_fig()`) |
| `lib_ml.R` | leakage-free pipeline, metrics, permutation |
| `01_prepare_data.R` | recoding from SPSS; sanity checks against the manuscript |
| `02_descriptives.R` | missingness; with vs without qEEG; cEEG-to-death timing; qEEG vs PMA, vigilance state and channels; epilepsy follow-up and attrition |
| `03_univariable.R` | Firth ORs per SD; BH-FDR; PMA- and vigilance-adjusted ORs |
| `04_legacy_bridge.R` | pipeline; same model with preprocessing inside folds; asymmetry artefact |
| `05_ml_main.R` | Clinical / EEG / combined; LightGBM and ridge; delta AUC; permutation tests |
| `06_ml_sensitivity.R` | S1-S8 sensitivity analyses; death by redirection of care; predictions vs channel count |
| `07_shap_stability.R` | SHAP rank stability across folds and bootstrap refits |
| `08_figures.R` | 
| `09_report.R` | 
| `run_all.R` | 

## Requirements

R ≥ 4.4 and the packages in `DESCRIPTION` (installed by
`R/install_requirements.R`). The exact versions behind the reported numbers are
in `output/tables/session_provenance.txt`. Optional: Ghostscript, which embeds
fonts in the PDF figures when cairo is unavailable.
