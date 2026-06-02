## ============================================================
## SENSITIVITY CHECK (Ben comment #1): Wider difficulty spread
##
## Re-runs ONLY the cross-population misalignment family with
## item difficulties drawn from a WIDER distribution:
##     b ~ N(0, 1.5^2)   instead of   b ~ N(0, 0.75^2)
##
## Tests the headline finding of Zhang, Liu, Molenaar & Domingue
## ("Realistic Simulation of Item Difficulties"): that difficulty
## SPREAD is the distributional feature that matters most.
##
## Question answered:
##   Under a wider difficulty spread, does
##     (a) scoring error still dominate (paper: 86-96%), and
##     (b) rescoring still drive CAT-R bias below 0.5% of MSE?
##
## SELF-CONTAINED. Does NOT source run_full_simulation_v2.R.
## Does NOT read or write output/data/res.rds. Nothing you
## already have is touched.
##
## The ONLY change from the main simulation is b_sd (0.75 -> 1.5).
## Everything else (100 items, 100 examinees, 50 reps, 20-item
## CAT, EAP N(0,1.5) prior, 20 rho levels) matches Study 1.
##
## Runtime: ~30-45 min (cross-population family only).
## USAGE: Rscript run_widespread_sensitivity.R
## ============================================================

library(tidyverse)

# ── Settings (match main sim EXCEPT b_sd) ─────────────────────
n_items     <- 100
n_examinees <- 100
n_reps      <- 50
cat_length  <- 20
prior_sd    <- 1.5
grid        <- seq(-6, 6, by = 0.05)

b_sd_wide   <- 1.5     # <<< THE ONLY CHANGE. Main sim uses 0.75.

rho_levels  <- seq(0.01, 0.99, length.out = 20)   # same grid as main sim

out_dir <- "output/data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Generation (b_sd is now a parameter) ──────────────────────
generate_pool <- function(n_items, b_sd) {
  tibble(
    item = 1:n_items,
    b = rnorm(n_items, 0, b_sd),
    a = pmax(0.5, exp(rnorm(n_items, 0.1, 0.3)))
  )
}
generate_examinees <- function(n) rnorm(n, 0, 1)

# ── IRT ───────────────────────────────────────────────────────
p_2pl <- function(theta, a, b) 1 / (1 + exp(-a * (theta - b)))
info_2pl <- function(theta, a, b) { p <- p_2pl(theta, a, b); a^2 * p * (1 - p) }
eap_score <- function(responses, a, b, grid, prior_sd) {
  ll <- sapply(grid, function(th) {
    p <- p_2pl(th, a, b); p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
    sum(responses * log(p) + (1 - responses) * log(1 - p))
  })
  lp <- ll + dnorm(grid, 0, prior_sd, log = TRUE)
  w  <- exp(lp - max(lp)); w <- w / sum(w)
  sum(grid * w)
}

# ── CAT + fixed scoring (identical to main sim) ───────────────
run_cat <- function(theta_true, pool, select_bank, score_bank,
                    cat_length, grid, prior_sd) {
  n_pool <- nrow(pool); used <- logical(n_pool)
  responses <- numeric(0); items_given <- integer(0); theta_hat <- 0
  for (k in seq_len(cat_length)) {
    info_vals <- info_2pl(theta_hat, select_bank$a, select_bank$b)
    info_vals[used] <- -Inf
    chosen <- which.max(info_vals); used[chosen] <- TRUE
    items_given <- c(items_given, chosen)
    p <- p_2pl(theta_true, pool$a[chosen], pool$b[chosen])
    responses <- c(responses, rbinom(1, 1, p))
    theta_hat <- eap_score(responses, select_bank$a[items_given],
                           select_bank$b[items_given], grid, prior_sd)
  }
  theta_final <- eap_score(responses, score_bank$a[items_given],
                           score_bank$b[items_given], grid, prior_sd)
  list(theta = theta_final, items = items_given, responses = responses)
}
run_fixed <- function(theta_true, pool, items_given, score_bank, grid, prior_sd) {
  responses <- rbinom(length(items_given), 1,
                      p_2pl(theta_true, pool$a[items_given], pool$b[items_given]))
  theta <- eap_score(responses, score_bank$a[items_given],
                     score_bank$b[items_given], grid, prior_sd)
  list(theta = theta)
}

# ── Cross-population misalignment (variance-preserving) ───────
# Identical to main sim's add_correlated(): mu and sd_b are read
# from the pool, so a wider pool flows through cleanly. No
# centering/standardising step silently undoes the spread.
add_correlated <- function(pool, rho) {
  pn <- pool; mu <- mean(pool$b); sd_b <- sd(pool$b)
  residual_sd <- sd_b * sqrt(max(1 - rho^2, 0))
  pn$b <- mu + rho * (pool$b - mu) + rnorm(nrow(pool), 0, residual_sd)
  pn
}

# ── Evaluate one cell (5 conditions; per-rep metrics) ─────────
evaluate_cell <- function(pool, misaligned_pool, n_examinees, n_reps,
                          cat_length, grid, prior_sd) {
  true_bank   <- list(a = pool$a, b = pool$b)
  noisy_bank  <- list(a = misaligned_pool$a, b = misaligned_pool$b)
  fixed_items <- order(pool$a, decreasing = TRUE)[1:cat_length]
  rows <- list(); idx <- 0
  for (rep in seq_len(n_reps)) {
    thetas <- generate_examinees(n_examinees)
    for (i in seq_along(thetas)) {
      th <- thetas[i]
      cat_o <- run_cat(th, pool, noisy_bank, noisy_bank, cat_length, grid, prior_sd)
      theta_r <- eap_score(cat_o$responses, true_bank$a[cat_o$items],
                           true_bank$b[cat_o$items], grid, prior_sd)
      cat_b <- run_cat(th, pool, true_bank, true_bank, cat_length, grid, prior_sd)
      fix   <- run_fixed(th, pool, fixed_items, true_bank, grid, prior_sd)
      full  <- run_fixed(th, pool, seq_len(n_items), true_bank, grid, prior_sd)
      base <- tibble(rep = rep, examinee = i, theta_true = th)
      idx <- idx+1; rows[[idx]] <- base %>% mutate(condition="CAT-O",    theta_hat=cat_o$theta)
      idx <- idx+1; rows[[idx]] <- base %>% mutate(condition="CAT-R",    theta_hat=theta_r)
      idx <- idx+1; rows[[idx]] <- base %>% mutate(condition="CAT-B",    theta_hat=cat_b$theta)
      idx <- idx+1; rows[[idx]] <- base %>% mutate(condition="Fixed-20", theta_hat=fix$theta)
      idx <- idx+1; rows[[idx]] <- base %>% mutate(condition="Full-100", theta_hat=full$theta)
    }
  }
  bind_rows(rows)
}

compute_metrics_by_rep <- function(cell_results) {
  cell_results %>%
    group_by(condition, rep) %>%
    summarise(
      bias = mean(theta_hat - theta_true),
      mse  = mean((theta_hat - theta_true)^2),
      .groups = "drop"
    )
}

# ══════════════════════════════════════════════════════════════
# RUN: cross-population family under WIDE spread
# ══════════════════════════════════════════════════════════════
set.seed(2026)
pool <- generate_pool(n_items, b_sd_wide)
cat("Wide-spread pool generated. b_sd target =", b_sd_wide, "\n")
cat("  observed SD(b) =", round(sd(pool$b), 3),
    "| range [", round(min(pool$b),2), ",", round(max(pool$b),2), "]\n\n")

set.seed(2027)
all_results <- list()
t0 <- Sys.time()
for (r in seq_along(rho_levels)) {
  rho <- rho_levels[r]
  cat(sprintf("  rho = %.3f (%d/%d) [%s]\n", rho, r, length(rho_levels),
              format(Sys.time() - t0)))
  mis  <- add_correlated(pool, rho)
  cell <- evaluate_cell(pool, mis, n_examinees, n_reps, cat_length, grid, prior_sd)
  metrics <- compute_metrics_by_rep(cell)
  all_results <- c(all_results, list(metrics %>% mutate(rho = rho)))
}
res <- bind_rows(all_results)

# ── Aggregate to cell level (same math as analysis script) ────
res_cell <- res %>%
  group_by(rho, condition) %>%
  summarise(
    cell_rmse     = sqrt(mean(mse)),
    cell_bias     = mean(bias),
    cell_mse      = mean(mse),
    cell_pct_bias = if_else(mean(mse) > 1e-10, mean(bias)^2 / mean(mse) * 100, 0),
    .groups = "drop"
  )

# ══════════════════════════════════════════════════════════════
# THE TWO NUMBERS BEN'S COMMENT NEEDS
# ══════════════════════════════════════════════════════════════
cat("\n\n========== WIDE-SPREAD SENSITIVITY RESULTS (b_sd =", b_sd_wide, ") ==========\n")

# (a) Scoring error as % of miscalibration-attributable error
scoring <- res_cell %>%
  select(rho, condition, cell_rmse) %>%
  pivot_wider(names_from = condition, values_from = cell_rmse) %>%
  mutate(
    total_error   = `CAT-O` - `CAT-B`,
    scoring_error = `CAT-O` - `CAT-R`,
    pct_scoring   = if_else(total_error > 1e-6, scoring_error / total_error * 100, NA_real_)
  )
mean_pct_scoring <- mean(scoring$pct_scoring, na.rm = TRUE)

cat("\n(a) Scoring error as % of miscalibration-attributable error:\n")
cat(sprintf("    Main sim (b_sd=0.75) reported: 87%% for cross-population.\n"))
cat(sprintf("    Wide spread (b_sd=1.5) mean : %.1f%%\n", mean_pct_scoring))

# (b) CAT-R bias as % of MSE
catr_bias <- res_cell %>% filter(condition == "CAT-R")
cat("\n(b) CAT-R bias as % of MSE:\n")
cat(sprintf("    Main sim reported: < 0.5%% (cross-pop max ~0.13%%).\n"))
cat(sprintf("    Wide spread mean : %.4f%%   max: %.4f%%\n",
            mean(catr_bias$cell_pct_bias), max(catr_bias$cell_pct_bias)))

# CAT-O bias for contrast
cato_bias <- res_cell %>% filter(condition == "CAT-O")
cat(sprintf("    (CAT-O bias for contrast: mean %.2f%%, max %.2f%%)\n",
            mean(cato_bias$cell_pct_bias), max(cato_bias$cell_pct_bias)))

cat("\n--- Verdict ---\n")
cat(sprintf("  Scoring dominates (>= ~85%%)?  %s\n",
            if_else(mean_pct_scoring >= 85, "YES", "CHECK")))
cat(sprintf("  CAT-R bias stays < 0.5%%?      %s\n",
            if_else(max(catr_bias$cell_pct_bias) < 0.5, "YES", "CHECK")))

# Save a small CSV so the numbers are recoverable later
write_csv(res_cell, file.path(out_dir, "sensitivity_widespread_cells.csv"))
write_csv(scoring,  file.path(out_dir, "sensitivity_widespread_scoring.csv"))
cat("\nSaved: output/data/sensitivity_widespread_cells.csv\n")
cat("       output/data/sensitivity_widespread_scoring.csv\n")
cat("\nDone!\n")