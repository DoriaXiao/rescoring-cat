## ============================================================
## Rerun stepwise simulation with more replications
##
## Purpose: Mike (Comment 10) flagged a "weird bump" on the left
## side of CAT-O in Figure 4. With 40 reps, CAT-O's step-2 RMSE
## (0.94) is non-monotonically higher than step-1 (0.93). All other
## conditions decline monotonically. We re-run with 200 reps to
## test whether the bump shrinks (small-sample artifact) or
## persists (structural EAP+miscalibration interaction at low n).
##
## All other settings match run_full_simulation_v2.R exactly:
##   - 2PL pool, 100 items, b ~ N(0, 0.75^2),
##     a ~ pmax(0.5, exp(rnorm(0.1, 0.3)))
##   - theta ~ N(0, 1), 300 examinees per rep
##   - rho = 0.75 (matches existing Figure 4)
##   - cat_length = 20 (record per-step RMSE)
##   - EAP scoring with N(0, 1.5) prior, grid [-6, 6] step 0.05
##   - Deterministic max-info selection
##
## Output: output/data/res_stepwise_v2.rds
##         columns: condition, step, rmse, bias, cor_r
##         (drop-in replacement for the original res_stepwise.rds;
##         identical schema. To use it in Figure 4, change line 20 of
##         run_simulation_analysis_v2.R from "res_stepwise.rds" to
##         "res_stepwise_v2.rds".)
##
## Runtime: ~1.5-2.5 hours (200 reps × 300 examinees × 80 EAP calls
##          per examinee = 4.8M EAP fits)
## ============================================================

library(tidyverse)
library(here)

# ── 0. Settings ───────────────────────────────────────────────

n_items      <- 100
n_examinees  <- 300            # match caption
n_reps       <- 200            # 5x more than original
cat_length   <- 20
prior_sd     <- 1.5
grid         <- seq(-6, 6, by = 0.05)
rho_target   <- 0.75           # match Figure 4

out_dir <- here("output", "data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ── 1. IRT helpers (copied from run_full_simulation_v2.R) ────

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


# ── 2. Stepwise CAT engine ────────────────────────────────────
# Returns per-step theta_hat under three scoring conditions:
# CAT-O: select & score with noisy bank
# CAT-R: select with noisy bank, score with true bank
# CAT-B: select & score with true bank (separate run, separate items)

run_cat_stepwise <- function(theta_true, true_bank, noisy_bank,
                             cat_length, grid, prior_sd) {
  n_pool      <- length(true_bank$b)

  # Two parallel CAT runs needed: one with noisy selection (gives O & R),
  # one with true selection (gives B).
  # Both use the same theta_true so they're paired.

  # ── Run 1: CAT with noisy selection bank ──
  used_n     <- logical(n_pool); items_n <- integer(0); resp_n <- numeric(0)
  th_o_step  <- numeric(cat_length)   # CAT-O step trajectory
  th_r_step  <- numeric(cat_length)   # CAT-R step trajectory
  th_hat_n   <- 0

  for (k in seq_len(cat_length)) {
    info_vals <- info_2pl(th_hat_n, noisy_bank$a, noisy_bank$b)
    info_vals[used_n] <- -Inf
    chosen <- which.max(info_vals)
    used_n[chosen] <- TRUE
    items_n <- c(items_n, chosen)

    p <- p_2pl(theta_true, true_bank$a[chosen], true_bank$b[chosen])
    resp_n <- c(resp_n, rbinom(1, 1, p))

    # Update for next-item selection (uses noisy bank)
    th_hat_n <- eap_score(resp_n, noisy_bank$a[items_n],
                          noisy_bank$b[items_n], grid, prior_sd)

    # Record CAT-O score (= the running noisy-bank EAP)
    th_o_step[k] <- th_hat_n

    # Record CAT-R score (rescore same responses with true bank)
    th_r_step[k] <- eap_score(resp_n, true_bank$a[items_n],
                              true_bank$b[items_n], grid, prior_sd)
  }

  # ── Run 2: CAT with true selection bank (CAT-B benchmark) ──
  used_b     <- logical(n_pool); items_b <- integer(0); resp_b <- numeric(0)
  th_b_step  <- numeric(cat_length)
  th_hat_b   <- 0

  for (k in seq_len(cat_length)) {
    info_vals <- info_2pl(th_hat_b, true_bank$a, true_bank$b)
    info_vals[used_b] <- -Inf
    chosen <- which.max(info_vals)
    used_b[chosen] <- TRUE
    items_b <- c(items_b, chosen)

    p <- p_2pl(theta_true, true_bank$a[chosen], true_bank$b[chosen])
    resp_b <- c(resp_b, rbinom(1, 1, p))

    th_hat_b <- eap_score(resp_b, true_bank$a[items_b],
                          true_bank$b[items_b], grid, prior_sd)
    th_b_step[k] <- th_hat_b
  }

  list(`CAT-O` = th_o_step, `CAT-R` = th_r_step, `CAT-B` = th_b_step)
}


# ── 3. Fixed-20 stepwise (separate; doesn't depend on selection) ──

run_fixed_stepwise <- function(theta_true, true_bank, fixed_items,
                               grid, prior_sd) {
  resp <- rbinom(length(fixed_items), 1,
                 p_2pl(theta_true, true_bank$a[fixed_items],
                       true_bank$b[fixed_items]))
  th_step <- numeric(length(fixed_items))
  for (k in seq_along(fixed_items)) {
    use <- fixed_items[1:k]
    th_step[k] <- eap_score(resp[1:k], true_bank$a[use],
                            true_bank$b[use], grid, prior_sd)
  }
  th_step
}


# ── 4. Run all reps ───────────────────────────────────────────

cat("=== Stepwise rerun: rho =", rho_target, "===\n")
cat("Reps:", n_reps, "| Examinees per rep:", n_examinees, "\n\n")

set.seed(2026)
shared_pool <- generate_pool(n_items)

# Generate one misaligned pool (rho=0.75) and use the same one across reps,
# matching how Figure 4 was originally produced.
set.seed(2027)
misaligned_pool <- add_correlated(shared_pool, rho_target)
true_bank  <- list(a = shared_pool$a,    b = shared_pool$b)
noisy_bank <- list(a = misaligned_pool$a, b = misaligned_pool$b)
fixed_items <- order(shared_pool$a, decreasing = TRUE)[1:cat_length]

# Storage: list of arrays (n_reps x cat_length) per condition
storage <- list(
  `CAT-O`    = matrix(NA_real_, n_reps * n_examinees, cat_length),
  `CAT-R`    = matrix(NA_real_, n_reps * n_examinees, cat_length),
  `CAT-B`    = matrix(NA_real_, n_reps * n_examinees, cat_length),
  `Fixed-20` = matrix(NA_real_, n_reps * n_examinees, cat_length)
)
theta_true_all <- numeric(n_reps * n_examinees)

t0 <- Sys.time()
row <- 0

for (r in seq_len(n_reps)) {
  if (r %% 10 == 0) {
    cat(sprintf("  rep %d/%d  [%s elapsed]\n",
                r, n_reps, format(Sys.time() - t0)))
  }

  set.seed(3000 + r)
  thetas <- rnorm(n_examinees)

  for (i in seq_along(thetas)) {
    row <- row + 1
    th <- thetas[i]
    theta_true_all[row] <- th

    cat_steps <- run_cat_stepwise(th, true_bank, noisy_bank,
                                  cat_length, grid, prior_sd)
    storage$`CAT-O`[row, ] <- cat_steps$`CAT-O`
    storage$`CAT-R`[row, ] <- cat_steps$`CAT-R`
    storage$`CAT-B`[row, ] <- cat_steps$`CAT-B`

    storage$`Fixed-20`[row, ] <- run_fixed_stepwise(th, true_bank,
                                                    fixed_items, grid, prior_sd)
  }
}

cat("\nTotal time:", format(Sys.time() - t0), "\n")


# ── 5. Compute per-step RMSE / bias / cor ─────────────────────

stepwise <- map_dfr(names(storage), function(cond) {
  mat <- storage[[cond]]
  tibble(
    condition = cond,
    step      = 1:cat_length,
    rmse      = sqrt(colMeans((mat - theta_true_all)^2)),
    bias      = colMeans(mat - theta_true_all),
    cor_r     = sapply(seq_len(cat_length), function(k) cor(mat[, k], theta_true_all))
  )
})

cat("\n=== Stepwise summary ===\n")
print(stepwise %>%
      select(condition, step, rmse) %>%
      pivot_wider(names_from = condition, values_from = rmse) %>%
      slice(1:5))

cat("\nDoes the bump persist?\n")
o_step1 <- stepwise$rmse[stepwise$condition == "CAT-O" & stepwise$step == 1]
o_step2 <- stepwise$rmse[stepwise$condition == "CAT-O" & stepwise$step == 2]
cat(sprintf("  CAT-O step 1: %.4f\n", o_step1))
cat(sprintf("  CAT-O step 2: %.4f\n", o_step2))
cat(sprintf("  Difference (step2 - step1): %+.4f\n", o_step2 - o_step1))
if (o_step2 > o_step1) {
  cat("  -> BUMP PERSISTS (structural). Add caption explanation.\n")
} else {
  cat("  -> BUMP GONE (was small-sample artifact).\n")
}


# ── 6. Save ───────────────────────────────────────────────────

saveRDS(stepwise, file.path(out_dir, "res_stepwise_v2.rds"))
cat("\nSaved:", file.path(out_dir, "res_stepwise_v2.rds"), "\n")
cat("(Update run_simulation_analysis_v2.R to read res_stepwise_v2.rds\n")
cat(" if you want to use this in the figure.)\n")
