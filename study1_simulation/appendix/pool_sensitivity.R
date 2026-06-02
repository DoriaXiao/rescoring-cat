## ============================================================
## Sensitivity Check: Pool-to-Pool Variation in Family 1
##
## Purpose: Demonstrate that the rescoring effect (CAT-O vs CAT-R)
## is robust across different item pools. Produces an appendix
## figure showing CAT-O and CAT-R mean curves with 10th-90th
## percentile bands across 10 independently generated pools.
##
## NOTE: This is a stripped-down Family 1 replication. We use
## n_reps = 20 (vs 50 in main sim) to keep runtime manageable;
## with 10 pools x 100 examinees x 20 reps = 20,000 scores per
## rho level per pool, which is more than sufficient for bands.
##
## Runtime: ~4-5 hours (roughly 5x original Family 1).
## Run independently of run_full_simulation_patched.R.
## ============================================================

library(tidyverse)

# ── Settings ──────────────────────────────────────────────────

n_items        <- 100
n_examinees    <- 100
n_reps         <- 20     # reduced from 50 for runtime
cat_length     <- 20
prior_sd       <- 1.5
grid           <- seq(-6, 6, by = 0.05)
n_pools        <- 10
rho_levels     <- seq(0.01, 0.99, length.out = 20)

out_dir <- "output/data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ── Helper functions (identical to main simulation) ────────────

generate_pool <- function(n_items) {
  tibble(
    item = 1:n_items,
    b = rnorm(n_items, 0, 0.75),
    a = pmax(0.5, exp(rnorm(n_items, 0.1, 0.3)))
  )
}

generate_examinees <- function(n) rnorm(n, 0, 1)

p_2pl <- function(theta, a, b) 1 / (1 + exp(-a * (theta - b)))

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

run_cat <- function(theta_true, pool, select_bank, score_bank,
                    cat_length, grid, prior_sd) {
  n_pool <- nrow(pool)
  used <- logical(n_pool)
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

add_correlated <- function(pool, rho) {
  pool_noisy <- pool
  mu <- mean(pool$b)
  sd_b <- sd(pool$b)
  residual_sd <- sd_b * sqrt(max(1 - rho^2, 0))
  pool_noisy$b <- mu + rho * (pool$b - mu) + rnorm(nrow(pool), 0, residual_sd)
  pool_noisy
}

evaluate_cell <- function(pool, misaligned_pool, n_examinees, n_reps,
                          cat_length, grid, prior_sd) {
  true_bank  <- list(a = pool$a, b = pool$b)
  noisy_bank <- list(a = misaligned_pool$a, b = misaligned_pool$b)

  results <- vector("list", n_reps * n_examinees * 3)
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

      base <- tibble(rep = rep, examinee = i, theta_true = th)
      idx <- idx + 1; results[[idx]] <- base %>% mutate(condition = "CAT-O", theta_hat = cat_o$theta)
      idx <- idx + 1; results[[idx]] <- base %>% mutate(condition = "CAT-R", theta_hat = theta_r)
      idx <- idx + 1; results[[idx]] <- base %>% mutate(condition = "CAT-B", theta_hat = cat_b$theta)
    }
  }

  bind_rows(results[1:idx]) %>%
    group_by(condition) %>%
    summarise(rmse = sqrt(mean((theta_hat - theta_true)^2)),
              .groups = "drop")
}


# ══════════════════════════════════════════════════════════════
# MAIN LOOP: 10 pools x 20 rho levels
# ══════════════════════════════════════════════════════════════

set.seed(3026)
all_results <- list()
t0 <- Sys.time()

for (p in seq_len(n_pools)) {
  cat(sprintf("=== Pool %d/%d [%s] ===\n", p, n_pools, format(Sys.time() - t0)))
  pool <- generate_pool(n_items)

  for (r in seq_along(rho_levels)) {
    rho <- rho_levels[r]
    if (r %% 5 == 0) cat(sprintf("  rho = %.3f (%d/%d)\n",
                                  rho, r, length(rho_levels)))

    misaligned <- add_correlated(pool, rho)
    metrics <- evaluate_cell(pool, misaligned, n_examinees, n_reps,
                             cat_length, grid, prior_sd)

    all_results <- c(all_results, list(
      metrics %>% mutate(pool_id = p, rho = rho)
    ))
  }
}

sens <- bind_rows(all_results)

cat("\nTotal time:", format(Sys.time() - t0), "\n")
cat("rows:", nrow(sens), "\n")

saveRDS(sens, file.path(out_dir, "res_sensitivity.rds"))
cat("Saved:", file.path(out_dir, "res_sensitivity.rds"), "\n")
