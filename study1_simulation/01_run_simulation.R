## ============================================================
## PATCHED CAT Rescoring Simulation (v2)
##
## Two fixes over the original run_full_simulation.R:
##   1. Shared item pool across all three families
##   2. Per-rep metric output (enables SE bands on figures)
##
## Output format: one row per (condition, rep) per cell.
## The analysis script aggregates to cell-level means and SEs.
##
## Does NOT re-run stepwise trajectory (already clean).
##
## USAGE
##   1. Back up: cp output/data/res.rds output/data/res.rds.backup
##   2. Run this script (~2-3 hours)
##   3. Run run_simulation_analysis_v2.R to regenerate figures
##
## Runtime: ~2-3 hours (families 1-3 only)
## ============================================================

library(tidyverse)

# ── 0. Settings (match paper exactly) ─────────────────────────

n_items      <- 100
n_examinees  <- 100
n_reps       <- 50
cat_length   <- 20
prior_sd     <- 1.5
grid         <- seq(-6, 6, by = 0.05)

out_dir <- "output/data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ── 1. Item & examinee generation ─────────────────────────────

generate_pool <- function(n_items) {
  tibble(
    item = 1:n_items,
    b = rnorm(n_items, 0, 0.75),
    a = pmax(0.5, exp(rnorm(n_items, 0.1, 0.3)))
  )
}

generate_examinees <- function(n) {
  rnorm(n, 0, 1)
}


# ── 2. IRT functions ──────────────────────────────────────────

p_2pl <- function(theta, a, b) {
  1 / (1 + exp(-a * (theta - b)))
}

info_2pl <- function(theta, a, b) {
  p <- p_2pl(theta, a, b)
  a^2 * p * (1 - p)
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


# ── 3. CAT engine ─────────────────────────────────────────────

run_cat <- function(theta_true, pool, select_bank, score_bank,
                    cat_length, grid, prior_sd) {
  n_pool    <- nrow(pool)
  used      <- logical(n_pool)
  responses <- numeric(0)
  items_given <- integer(0)
  theta_hat <- 0

  for (k in seq_len(cat_length)) {
    info_vals <- info_2pl(theta_hat, select_bank$a, select_bank$b)
    info_vals[used] <- -Inf
    chosen <- which.max(info_vals)
    used[chosen] <- TRUE
    items_given <- c(items_given, chosen)

    p <- p_2pl(theta_true, pool$a[chosen], pool$b[chosen])
    resp <- rbinom(1, 1, p)
    responses <- c(responses, resp)

    theta_hat <- eap_score(responses, select_bank$a[items_given],
                           select_bank$b[items_given], grid, prior_sd)
  }

  theta_final <- eap_score(responses, score_bank$a[items_given],
                           score_bank$b[items_given], grid, prior_sd)

  list(theta = theta_final, items = items_given, responses = responses)
}

run_fixed <- function(theta_true, pool, items_given, score_bank,
                      grid, prior_sd) {
  responses <- rbinom(length(items_given), 1,
                      p_2pl(theta_true, pool$a[items_given], pool$b[items_given]))
  theta <- eap_score(responses, score_bank$a[items_given],
                     score_bank$b[items_given], grid, prior_sd)
  list(theta = theta, items = items_given, responses = responses)
}


# ── 4. Misalignment functions ─────────────────────────────────

add_noise <- function(pool, sigma) {
  pool_noisy <- pool
  if (sigma == 0) return(pool_noisy)
  pool_noisy$b <- pool$b + rnorm(nrow(pool), 0, sigma)
  pool_noisy
}

add_shift_noise <- function(pool, delta, sigma) {
  pool_noisy <- pool
  pool_noisy$b <- pool$b + delta + rnorm(nrow(pool), 0, sigma)
  pool_noisy
}

add_correlated <- function(pool, rho) {
  pool_noisy <- pool
  mu <- mean(pool$b)
  sd_b <- sd(pool$b)
  residual_sd <- sd_b * sqrt(max(1 - rho^2, 0))
  pool_noisy$b <- mu + rho * (pool$b - mu) + rnorm(nrow(pool), 0, residual_sd)
  pool_noisy
}


# ── 5. Evaluate one cell (unchanged from original) ────────────

evaluate_cell <- function(pool, misaligned_pool, n_examinees, n_reps,
                          cat_length, grid, prior_sd) {
  true_bank  <- list(a = pool$a, b = pool$b)
  noisy_bank <- list(a = misaligned_pool$a, b = misaligned_pool$b)
  fixed_items <- order(pool$a, decreasing = TRUE)[1:cat_length]

  results <- vector("list", n_reps * n_examinees * 5)
  idx <- 0

  for (rep in seq_len(n_reps)) {
    thetas <- generate_examinees(n_examinees)

    for (i in seq_along(thetas)) {
      th <- thetas[i]

      cat_o <- run_cat(th, pool, noisy_bank, noisy_bank,
                       cat_length, grid, prior_sd)

      theta_r <- eap_score(cat_o$responses, true_bank$a[cat_o$items],
                           true_bank$b[cat_o$items], grid, prior_sd)

      cat_b <- run_cat(th, pool, true_bank, true_bank,
                       cat_length, grid, prior_sd)

      fix <- run_fixed(th, pool, fixed_items, true_bank, grid, prior_sd)

      full <- run_fixed(th, pool, seq_len(n_items), true_bank, grid, prior_sd)

      base <- tibble(rep = rep, examinee = i, theta_true = th)
      idx <- idx + 1; results[[idx]] <- base %>% mutate(condition = "CAT-O", theta_hat = cat_o$theta)
      idx <- idx + 1; results[[idx]] <- base %>% mutate(condition = "CAT-R", theta_hat = theta_r)
      idx <- idx + 1; results[[idx]] <- base %>% mutate(condition = "CAT-B", theta_hat = cat_b$theta)
      idx <- idx + 1; results[[idx]] <- base %>% mutate(condition = "Fixed-20", theta_hat = fix$theta)
      idx <- idx + 1; results[[idx]] <- base %>% mutate(condition = "Full-100", theta_hat = full$theta)
    }
  }

  bind_rows(results[1:idx])
}


# ── 6. Compute metrics PER REP (the key change from v1) ──────
#
# Returns one row per (condition, rep), preserving replication
# structure for SE band computation downstream.
# Cell-level aggregates (for paper numbers) are recovered in
# the analysis script by averaging across reps.

compute_metrics_by_rep <- function(cell_results) {
  cell_results %>%
    group_by(condition, rep) %>%
    summarise(
      rmse = sqrt(mean((theta_hat - theta_true)^2)),
      bias = mean(theta_hat - theta_true),
      cor_r = cor(theta_hat, theta_true),
      mse  = mean((theta_hat - theta_true)^2),
      n    = n(),
      .groups = "drop"
    )
}


# ══════════════════════════════════════════════════════════════
# 7. GENERATE SHARED POOL
# ══════════════════════════════════════════════════════════════

set.seed(2026)
shared_pool <- generate_pool(n_items)

cat("Shared pool generated.\n")
cat("  b range: [", round(min(shared_pool$b), 2), ",",
    round(max(shared_pool$b), 2), "]\n")
cat("  a range: [", round(min(shared_pool$a), 2), ",",
    round(max(shared_pool$a), 2), "]\n\n")


# ══════════════════════════════════════════════════════════════
# 8. RUN ALL THREE FAMILIES
# ══════════════════════════════════════════════════════════════

all_results <- list()
t0 <- Sys.time()

# --- Family 1: Cross-population misalignment (rho) ---
set.seed(2027)
cat("=== Family 1: Cross-population misalignment ===\n")
rho_levels <- seq(0.01, 0.99, length.out = 20)

for (r in seq_along(rho_levels)) {
  rho <- rho_levels[r]
  cat(sprintf("  rho = %.3f (%d/%d) [%s]\n", rho, r, length(rho_levels),
              format(Sys.time() - t0)))

  misaligned <- add_correlated(shared_pool, rho)

  cell <- evaluate_cell(shared_pool, misaligned, n_examinees, n_reps,
                        cat_length, grid, prior_sd)
  metrics <- compute_metrics_by_rep(cell)

  all_results <- c(all_results, list(
    metrics %>% mutate(mode = "correlated", rho = rho,
                       sigma = NA_real_, delta = NA_real_)
  ))
}


# --- Family 2: Random noise (sigma) ---
set.seed(2028)
cat("\n=== Family 2: Random noise ===\n")
sigma_levels <- c(0, seq(0.05, 1.0, length.out = 19))

for (s in seq_along(sigma_levels)) {
  sig <- sigma_levels[s]
  cat(sprintf("  sigma = %.3f (%d/%d) [%s]\n", sig, s, length(sigma_levels),
              format(Sys.time() - t0)))

  misaligned <- add_noise(shared_pool, sig)

  cell <- evaluate_cell(shared_pool, misaligned, n_examinees, n_reps,
                        cat_length, grid, prior_sd)
  metrics <- compute_metrics_by_rep(cell)

  all_results <- c(all_results, list(
    metrics %>% mutate(mode = "random_noise", rho = NA_real_,
                       sigma = sig, delta = NA_real_)
  ))
}


# --- Family 3: Structured shift + noise ---
set.seed(2029)
cat("\n=== Family 3: Structured shifts ===\n") 
delta_levels <- c(-0.3, 0.3, 1.0)
sigma_shift_levels <- seq(0.05, 1.0, length.out = 10)

for (d in delta_levels) {
  for (s in seq_along(sigma_shift_levels)) {
    sig <- sigma_shift_levels[s]
    cat(sprintf("  delta = %.1f, sigma = %.3f (%d/%d) [%s]\n",
                d, sig, s, length(sigma_shift_levels),
                format(Sys.time() - t0)))

    misaligned <- add_shift_noise(shared_pool, d, sig)

    cell <- evaluate_cell(shared_pool, misaligned, n_examinees, n_reps,
                          cat_length, grid, prior_sd)
    metrics <- compute_metrics_by_rep(cell)

    all_results <- c(all_results, list(
      metrics %>% mutate(mode = "structured_shift", rho = NA_real_,
                         sigma = sig, delta = d)
    ))
  }
}


# ══════════════════════════════════════════════════════════════
# 9. SAVE
# ══════════════════════════════════════════════════════════════

res <- bind_rows(all_results)

cat("\n\n========== SUMMARY ==========\n")
cat("Total time:", format(Sys.time() - t0), "\n")
cat("res.rds:", nrow(res), "rows\n")
cat("  (", n_distinct(res$rep), "reps per cell,",
    n_distinct(paste(res$mode, res$rho, res$sigma, res$delta)), "cells )\n")
count(res, mode) %>% print()

saveRDS(res, file.path(out_dir, "res.rds"))
cat("\nSaved:", file.path(out_dir, "res.rds"), "\n")
cat("(res_stepwise.rds unchanged)\n")


# ══════════════════════════════════════════════════════════════
# 10. QUICK CHECKS (aggregated from per-rep data)
# ══════════════════════════════════════════════════════════════

# Recover cell-level aggregates (same math as old compute_metrics)
res_cell <- res %>%
  group_by(mode, rho, sigma, delta, condition) %>%
  summarise(
    rmse     = sqrt(mean(mse)),
    bias     = mean(bias),
    mse      = mean(mse),
    bias_sq  = mean(bias)^2,
    pct_bias = if_else(mean(mse) > 1e-10, mean(bias)^2 / mean(mse) * 100, 0),
    .groups  = "drop"
  )

cat("\n-- Scoring error as % of total error, by family --\n")
res_cell %>%
  select(mode, rho, sigma, delta, condition, rmse) %>%
  pivot_wider(names_from = condition, values_from = rmse) %>%
  mutate(
    total_error   = `CAT-O` - `CAT-B`,
    scoring_error = `CAT-O` - `CAT-R`,
    pct_scoring   = if_else(total_error > 1e-6,
                            scoring_error / total_error * 100, NA_real_)
  ) %>%
  group_by(mode) %>%
  summarise(mean_pct_scoring = mean(pct_scoring, na.rm = TRUE),
            .groups = "drop") %>%
  print()

cat("\n-- Bias decomposition (Table 2 check) --\n")
res_cell %>%
  filter(condition %in% c("CAT-O", "CAT-R")) %>%
  group_by(mode, condition) %>%
  summarise(
    mean_pct_bias = mean(pct_bias),
    max_pct_bias  = max(pct_bias),
    .groups = "drop"
  ) %>%
  print()

cat("\n-- Max CAT-R pct_bias across all families --\n")
cat(sprintf("  %.4f%%\n",
            max(res_cell$pct_bias[res_cell$condition == "CAT-R"])))

cat("\nDone!\n")
