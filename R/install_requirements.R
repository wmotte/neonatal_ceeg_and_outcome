#!/usr/bin/env Rscript
#
# Installs the R packages needed by R/run_all.R (see DESCRIPTION).
# Exact versions used for the reported results: output/tables/session_provenance.txt
# and renv.lock. To restore those exact versions instead: renv::restore().
#
# Willem M. (Wim) Otte, w.m.otte@umcutrecht.nl
# ==============================================================================
pkgs <- c("tidyverse", "haven", "lightgbm", "glmnet", "logistf", "future", "furrr",
          "scales", "ragg", "writexl")
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, requireNamespace))
message("All required packages available: ", paste(pkgs, collapse = ", "))
# Optional system tool: Ghostscript (embeds fonts in PDF figures when cairo is
# unavailable). macOS: brew install ghostscript
