## ============================================================
## Paper Analysis 1: Bootstrap CIs on Cross-Task Correlation Diffs
## Paper Analysis 2: Bias Decomposition from Simulation
##
## Run from: ~/Documents/levante-runs/levante-runs/rescoring_simulation
## Requires: output/data/comparison.csv, output/data/res.rds
## ============================================================

library(tidyverse)

out_dir <- "output/tables"

# ══════════════════════════════════════════════════════════════
# ANALYSIS 1: Bootstrap CIs on cross-task correlation changes
# ══════════════════════════════════════════════════════════════

cat("=== Bootstrap CIs on cross-task validity changes ===\n\n")

comparison <- read_csv("output/data/comparison.csv", show_col_types = FALSE)

wide <- comparison %>%
  select(user_id, task_id, theta_old, theta_new_1pl) %>%
  pivot_wider(names_from = task_id,
              values_from = c(theta_old, theta_new_1pl))

task_pairs <- list(
  c("egma-math", "matrix-reasoning"),
  c("egma-math", "mental-rotation"),
  c("matrix-reasoning", "mental-rotation"),
  c("egma-math", "same-different-selection"),
  c("matrix-reasoning", "same-different-selection")
)

set.seed(2026)
n_boot <- 2000

boot_results <- map_dfr(task_pairs, function(pair) {
  o1 <- paste0("theta_old_", pair[1])
  o2 <- paste0("theta_old_", pair[2])
  n1 <- paste0("theta_new_1pl_", pair[1])
  n2 <- paste0("theta_new_1pl_", pair[2])

  d <- wide %>%
    filter(!is.na(.data[[o1]]), !is.na(.data[[o2]]),
           !is.na(.data[[n1]]), !is.na(.data[[n2]]))

  n <- nrow(d)
  r_old_obs  <- cor(d[[o1]], d[[o2]])
  r_new_obs  <- cor(d[[n1]], d[[n2]])
  delta_obs  <- r_new_obs - r_old_obs

  # Bootstrap
  deltas <- replicate(n_boot, {
    idx <- sample(n, n, replace = TRUE)
    r_old_b <- cor(d[[o1]][idx], d[[o2]][idx])
    r_new_b <- cor(d[[n1]][idx], d[[n2]][idx])
    r_new_b - r_old_b
  })

  tibble(
    pair      = paste0(pair[1], " x ", pair[2]),
    n         = n,
    r_old     = r_old_obs,
    r_new     = r_new_obs,
    delta     = delta_obs,
    ci_lower  = quantile(deltas, 0.025),
    ci_upper  = quantile(deltas, 0.975),
    p_value   = mean(deltas <= 0) * 2  # two-sided
  )
})

cat("Results (2000 bootstrap resamples):\n")
print(boot_results, width = 120)

# Interpretation helper
boot_results <- boot_results %>%
  mutate(significant = ci_lower > 0 | ci_upper < 0)

cat("\nAll significant?\n")
print(boot_results %>% select(pair, delta, ci_lower, ci_upper, significant))

write_csv(boot_results, file.path(out_dir, "results_cross_task_bootstrap.csv"))
cat("\nSaved results_cross_task_bootstrap.csv\n")


# ══════════════════════════════════════════════════════════════
# ANALYSIS 2: Bias decomposition from simulation
# ══════════════════════════════════════════════════════════════

cat("\n\n=== Simulation Bias Decomposition ===\n\n")

res <- readRDS("output/data/res.rds")

# Check that bias column exists
stopifnot("bias" %in% names(res))

# Create condition labels (same as plot_metric_accuracy_v4)
label_design <- function(design, select_bank, score_with) {
  dplyr::case_when(
    design == "CAT"   & select_bank == "bhat_cat"  & score_with == "noisy" ~ "CAT-O",
    design == "CAT"   & select_bank == "bhat_cat"  & score_with == "true"  ~ "CAT-R",
    design == "CAT"   & select_bank == "bhat_true" & score_with == "true"  ~ "CAT-B",
    design == "Fixed" & score_with == "true" ~ "Fixed-20",
    design == "Full"  & score_with == "true" ~ "Full-100",
    TRUE ~ NA_character_
  )
}

# --- Correlated family: bias and RMSE by rho ---
bias_rho <- res %>%
  filter(mode == "correlated") %>%
  mutate(condition = label_design(design, select_bank, score_with)) %>%
  filter(!is.na(condition)) %>%
  group_by(condition, rho) %>%
  summarise(
    mean_bias = mean(bias),
    sd_bias   = sd(bias),
    mean_rmse = mean(rmse),
    # Decompose: RMSE^2 = bias^2 + variance
    bias_sq   = mean(bias)^2,
    variance  = mean(rmse)^2 - mean(bias)^2,
    pct_bias  = mean(bias)^2 / mean(rmse)^2 * 100,  # % of MSE due to bias
    n_reps    = n(),
    .groups   = "drop"
  )

cat("-- Correlated family: Bias by rho (selected rho values) --\n")
bias_rho %>%
  filter(rho %in% c(0, 0.25, 0.5, 0.75, 0.9)) %>%
  select(condition, rho, mean_bias, mean_rmse, pct_bias) %>%
  arrange(rho, condition) %>%
  print(n = 30)

write_csv(bias_rho, file.path(out_dir, "results_sim_bias_rho.csv"))

# --- Additive family: bias and RMSE by sigma ---
bias_sigma <- res %>%
  filter(mode == "additive") %>%
  mutate(condition = label_design(design, select_bank, score_with)) %>%
  filter(!is.na(condition)) %>%
  group_by(condition, sigma, bias_type, b_bias) %>%
  summarise(
    mean_bias = mean(bias),
    mean_rmse = mean(rmse),
    pct_bias  = mean(bias)^2 / mean(rmse)^2 * 100,
    n_reps    = n(),
    .groups   = "drop"
  )

cat("\n-- Additive family: Bias by sigma (random noise, selected sigma) --\n")
bias_sigma %>%
  filter(bias_type == "random", sigma %in% c(0, 0.25, 0.5, 0.75, 1.0)) %>%
  select(condition, sigma, mean_bias, mean_rmse, pct_bias) %>%
  arrange(sigma, condition) %>%
  print(n = 30)

write_csv(bias_sigma, file.path(out_dir, "results_sim_bias_sigma.csv"))

# --- Summary for paper text ---
cat("\n-- Paper summary: CAT-O vs CAT-R bias at low rho --\n")
bias_rho %>%
  filter(condition %in% c("CAT-O", "CAT-R"), rho < 0.15) %>%
  select(condition, rho, mean_bias, mean_rmse, pct_bias) %>%
  print()

cat("\n-- Paper summary: how much of CAT-O's RMSE is bias vs variance? --\n")
bias_rho %>%
  filter(condition == "CAT-O") %>%
  summarise(
    avg_pct_bias_low_rho  = mean(pct_bias[rho < 0.3]),
    avg_pct_bias_mid_rho  = mean(pct_bias[rho >= 0.3 & rho < 0.7]),
    avg_pct_bias_high_rho = mean(pct_bias[rho >= 0.7])
  ) %>%
  print()

cat("\nDone! Saved:\n")
cat("  output/tables/results_cross_task_bootstrap.csv\n")
cat("  output/tables/results_sim_bias_rho.csv\n")
cat("  output/tables/results_sim_bias_sigma.csv\n")

