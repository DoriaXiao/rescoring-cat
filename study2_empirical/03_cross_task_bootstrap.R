## ============================================================
## Bootstrap CIs on Cross-Task Convergent-Validity Changes
##
## Bootstraps 95% CIs for the change in cross-task correlation
## (new vs old parameters) reported in the validation section.
##
## Run from the repo root (open rescoring-cat.Rproj).
## Requires: output/data/comparison.csv  (from 01_rescoring_analysis.R;
##             child-level, NOT bundled in the repo)
##
## NOTE: the simulation bias decomposition (Table 2 / Sec 2.4.3) is
## produced by study1_simulation/04_analysis_figures.R, which reads the
## current per-rep res.rds schema. The earlier combined script also did
## a bias decomposition keyed to the *old* res.rds schema; that stale
## half has been removed here and is preserved in archive/run_paper_analyses.R.
## ============================================================

library(tidyverse)
library(here)   # paths resolve relative to repo root (.here marker)

out_dir <- here("output", "tables")

# ══════════════════════════════════════════════════════════════
# ANALYSIS 1: Bootstrap CIs on cross-task correlation changes
# ══════════════════════════════════════════════════════════════

cat("=== Bootstrap CIs on cross-task validity changes ===\n\n")

comparison <- read_csv(here("output", "data", "comparison.csv"), show_col_types = FALSE)

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

cat("\nDone! Saved:\n")
cat("  output/tables/results_cross_task_bootstrap.csv\n")
