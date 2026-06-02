## ============================================================
## Run the full simulation and save res.rds
## Takes ~30 min with n_replicates=50, n_persons=100
## ============================================================

library(dplyr); library(ggplot2); library(purrr); library(tidyr)
setwd("~/Documents/levante-runs/levante-runs/rescoring_simulation")
# Source core functions (do NOT source plot_traj.R here — it overrides
# add_calibration_error with err_mode, which breaks run_combined_simulation)
source("R/simulation_source.R")

set.seed(2026)

res <- run_combined_simulation(
  sigma_draw  = "range", sigma_range = c(0, 1), n_sigma = 20,
  rho_draw    = "range", rho_range   = c(0, 1), n_rho = 20,
  b_sd_vec    = 0.75,
  stopping_rules = c("SEM", "fixed", "full"),
  b_bias_vec  = c(-0.3, 0.3, 1.0),
  n_persons   = 100,
  n_replicates = 50,
  enforce_max_length = TRUE,
  select_bank_opts = c("bhat_cat", "bhat_true"),
  modes       = c("correlated", "additive"),
  standardize_bprime = TRUE,
  N_cat = 150, N_rescore = 300,
  ignore_sigma_in_correlated = TRUE,
  sigma_anchor_for_correlated = 0
)

saveRDS(res, file = "res.rds")
cat("Saved res.rds:", nrow(res), "rows\n")

# Quick sanity checks
count(res, mode, bias_type, stopping_rule, select_bank)
