## ============================================================
## Perturbed-parameter sensitivity simulation
##
## Purpose: Address Mike's Slack ask #3 / Claude-as-reviewer's
## critique that the main simulation gives rescoring an unfairly
## clean win because CAT-R uses the TRUE difficulty parameters.
##
## In practice, "improved" rescoring parameters are estimates with
## their own residual error. This script adds a CAT-Rhat condition
## where rescoring uses bhat = b + eta, with
##     eta ~ N(0, sd(b)/sqrt(N_rescore))
## i.e., the standard noise scale from a finite calibration sample.
##
## We sweep three N_rescore values to show how rescoring quality
## degrades with smaller calibration samples:
##   N_rescore = 100, 300, 1000
## (300 matches the default in estimate_banks() in simulation_source.R)
##
## Restricted to cross-population misalignment (Mike: most informative
## family; the structured-shift bias claim doesn't depend on this
## because shift is fully identifiable from the offset).
##
## All other settings match run_full_simulation_v2.R exactly:
##   - 2PL pool, 100 items, b ~ N(0, 0.75^2),
##     a ~ pmax(0.5, exp(rnorm(0.1, 0.3)))
##   - theta ~ N(0, 1), 100 examinees per rep, 50 reps
##   - cat_length = 20
##   - EAP scoring with N(0, 1.5) prior, grid [-6, 6] step 0.05
##   - SHARED pool across all rho/N cells (matches main sim)
##
## Three rho values: 0.34 (matches LEVANTE matrix reasoning),
##                   0.50 (intermediate),
##                   0.75 (matches Figure 4 stepwise; LEVANTE mid)
##
## Output: output/data/res_perturbed.rds
##         columns: condition, rep, examinee, theta_true, theta_hat,
##                  rho, N_rescore
##
## Runtime estimate: ~30-65 minutes
##   (3 rho × 3 N × 50 reps × 100 examinees ~ 45k examinee-runs;
##    each runs 2 CATs (CAT-O, CAT-B) of 20 items + 2 final EAPs
##    = ~44 EAP calls per examinee = ~2M EAP calls total)
## ============================================================

library(tidyverse)
library(here)

# ── 0. Settings (match main sim) ──────────────────────────────

n_items      <- 100
n_examinees  <- 100
n_reps       <- 50
cat_length   <- 20
prior_sd     <- 1.5
grid         <- seq(-6, 6, by = 0.05)

rho_levels   <- c(0.34, 0.50, 0.75)
N_rescore_levels <- c(100, 300, 1000)

out_dir <- here("output", "data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ── 1. IRT helpers (copied from run_full_simulation_v2.R) ─────

p_2pl <- function(theta, a, b) 1 / (1 + exp(-a * (theta - b)))

info_2pl <- function(theta, a, b) {
  p <- p_2pl(theta, a, b); a^2 * p * (1 - p)
}

eap_score <- function(responses, a, b, grid, prior_sd) {
  ll <- sapply(grid, function(th) {
    p <- p_2pl(th, a, b)
    p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
    sum(responses * log(p) + (1 - responses) * log(1 - p))
  })
  lp <- ll + dnorm(grid, 0, prior_sd, log = TRUE)
  w  <- exp(lp - max(lp)); w <- w / sum(w)
  sum(grid * w)
}

generate_pool <- function(n_items) {
  tibble(
    item = 1:n_items,
    b = rnorm(n_items, 0, 0.75),
    a = pmax(0.5, exp(rnorm(n_items, 0.1, 0.3)))
  )
}

add_correlated <- function(pool, rho) {
  pool_noisy <- pool
  mu <- mean(pool$b); sd_b <- sd(pool$b)
  residual_sd <- sd_b * sqrt(max(1 - rho^2, 0))
  pool_noisy$b <- mu + rho * (pool$b - mu) + rnorm(nrow(pool), 0, residual_sd)
  pool_noisy
}

# Estimate b from a finite-N calibration: matches simulation_source.R
# `estimate_bank` helper:
#     bhat = b + N(0, sd(b)/sqrt(N))
estimate_b <- function(b_true, N) {
  se <- sd(b_true) / sqrt(max(1, N))
  b_true + rnorm(length(b_true), 0, se)
}


# ── 2. CAT engine (matches run_full_simulation_v2.R) ──────────

run_cat <- function(theta_true, pool, select_bank, score_bank,
                    cat_length, grid, prior_sd) {
  n_pool <- nrow(pool)
  used <- logical(n_pool); responses <- numeric(0); items_given <- integer(0)
  theta_hat <- 0

  for (k in seq_len(cat_length)) {
    info_vals <- info_2pl(theta_hat, select_bank$a, select_bank$b)
    info_vals[used] <- -Inf
    chosen <- which.max(info_vals)
    used[chosen] <- TRUE
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


# ── 3. Evaluate one cell (rho x N_rescore) ────────────────────

evaluate_cell <- function(pool, misaligned_pool, N_rescore,
                          n_examinees, n_reps, cat_length,
                          grid, prior_sd) {
  true_bank  <- list(a = pool$a, b = pool$b)
  noisy_bank <- list(a = misaligned_pool$a, b = misaligned_pool$b)

  # Cell-level estimated rescoring bank: bh = b + eta, drawn ONCE per cell.
  # This is the correct SYMMETRIC treatment relative to the main simulation,
  # which holds the misaligned pool fixed across all reps within a (rho) cell.
  # Both the operational miscalibration (b' = misaligned pool) and the
  # rescoring imperfection (bh = b + eta) are realized once per cell, so
  # all reps in the cell share both. CAT-O variation comes from theta and
  # response draws; CAT-Rhat variation does too, with the same rescoring bh
  # applied throughout — matching how a real recalibration would be deployed.
  bh <- estimate_b(true_bank$b, N_rescore)
  rescore_bank_est <- list(a = true_bank$a, b = bh)  # a kept exact
                                                     # (matches estimate_banks
                                                     # in simulation_source.R)

  results <- vector("list", n_reps * n_examinees * 4)
  idx <- 0

  for (rep in seq_len(n_reps)) {
    set.seed(5000 + rep + as.integer(N_rescore))  # per-rep reproducibility

    thetas <- rnorm(n_examinees)

    for (i in seq_along(thetas)) {
      th <- thetas[i]

      # CAT-O: select & score with noisy bank (operational)
      cat_o <- run_cat(th, pool, noisy_bank, noisy_bank,
                       cat_length, grid, prior_sd)

      # CAT-R (TRUE rescore): score with TRUE bank — current paper's claim
      theta_r_true <- eap_score(cat_o$responses,
                                true_bank$a[cat_o$items],
                                true_bank$b[cat_o$items],
                                grid, prior_sd)

      # CAT-Rhat (PERTURBED rescore): score with bh = b + eta
      # This is what's actually feasible in practice
      theta_r_hat <- eap_score(cat_o$responses,
                               rescore_bank_est$a[cat_o$items],
                               rescore_bank_est$b[cat_o$items],
                               grid, prior_sd)

      # CAT-B: select & score with true bank (counterfactual benchmark)
      cat_b <- run_cat(th, pool, true_bank, true_bank,
                       cat_length, grid, prior_sd)

      base <- tibble(rep = rep, examinee = i, theta_true = th,
                     N_rescore = N_rescore)
      idx <- idx + 1; results[[idx]] <- base %>% mutate(condition = "CAT-O",      theta_hat = cat_o$theta)
      idx <- idx + 1; results[[idx]] <- base %>% mutate(condition = "CAT-R",      theta_hat = theta_r_true)
      idx <- idx + 1; results[[idx]] <- base %>% mutate(condition = "CAT-Rhat",   theta_hat = theta_r_hat)
      idx <- idx + 1; results[[idx]] <- base %>% mutate(condition = "CAT-B",      theta_hat = cat_b$theta)
    }
  }

  bind_rows(results[1:idx])
}


# ── 4. Run all cells ──────────────────────────────────────────

cat("=== Perturbed-parameter sensitivity sim ===\n")
cat("rho levels:", paste(rho_levels, collapse = ", "), "\n")
cat("N_rescore levels:", paste(N_rescore_levels, collapse = ", "), "\n")
cat("Reps per cell:", n_reps, "| Examinees per rep:", n_examinees, "\n\n")

set.seed(2026)
shared_pool <- generate_pool(n_items)

all_results <- list()
t0 <- Sys.time()
total_cells <- length(rho_levels) * length(N_rescore_levels)
cell_counter <- 0

for (rho in rho_levels) {
  set.seed(2027 + as.integer(rho * 1000))
  misaligned_pool <- add_correlated(shared_pool, rho)

  for (N_r in N_rescore_levels) {
    cell_counter <- cell_counter + 1
    cat(sprintf("[%d/%d] rho = %.2f, N_rescore = %d  [%s]\n",
                cell_counter, total_cells, rho, N_r,
                format(Sys.time() - t0)))

    cell <- evaluate_cell(shared_pool, misaligned_pool, N_r,
                          n_examinees, n_reps, cat_length,
                          grid, prior_sd)
    cell$rho <- rho
    all_results <- c(all_results, list(cell))
  }
}

res_perturbed <- bind_rows(all_results)
cat("\nTotal time:", format(Sys.time() - t0), "\n")
cat("Rows:", nrow(res_perturbed), "\n")


# ── 5. Quick summary ──────────────────────────────────────────

summ <- res_perturbed %>%
  group_by(rho, N_rescore, condition) %>%
  summarise(
    rmse = sqrt(mean((theta_hat - theta_true)^2)),
    bias = mean(theta_hat - theta_true),
    bias_pct_mse = mean(theta_hat - theta_true)^2 /
                   mean((theta_hat - theta_true)^2) * 100,
    .groups = "drop"
  )

cat("\n=== Summary: RMSE by rho x N_rescore x condition ===\n")
summ %>%
  select(rho, N_rescore, condition, rmse) %>%
  pivot_wider(names_from = condition, values_from = rmse) %>%
  arrange(rho, N_rescore) %>%
  print(n = 100)

cat("\n=== Summary: Bias as %% of MSE ===\n")
summ %>%
  select(rho, N_rescore, condition, bias_pct_mse) %>%
  pivot_wider(names_from = condition, values_from = bias_pct_mse) %>%
  arrange(rho, N_rescore) %>%
  print(n = 100)

cat("\n=== Key claim check ===\n")
cat("Is CAT-Rhat consistently better than CAT-O across all cells?\n")
check <- summ %>%
  select(rho, N_rescore, condition, rmse) %>%
  pivot_wider(names_from = condition, values_from = rmse) %>%
  mutate(rhat_beats_o = `CAT-Rhat` < `CAT-O`,
         rhat_close_to_r = abs(`CAT-Rhat` - `CAT-R`) < abs(`CAT-Rhat` - `CAT-O`))
print(check %>% select(rho, N_rescore, rhat_beats_o, rhat_close_to_r),
      n = 100)

# ── 6. Save ───────────────────────────────────────────────────

saveRDS(res_perturbed, file.path(out_dir, "res_perturbed.rds"))
saveRDS(summ,          file.path(out_dir, "res_perturbed_summary.rds"))
cat("\nSaved:\n  ", file.path(out_dir, "res_perturbed.rds"),
    "\n  ", file.path(out_dir, "res_perturbed_summary.rds"), "\n")

cat("\nNext step: run perturbed-figure helper (or eyeball summary) and\n")
cat("decide whether to (a) keep perturbed sim as appendix-only sensitivity\n")
cat("check, or (b) replace CAT-R in main figures with CAT-Rhat at N=300.\n")
