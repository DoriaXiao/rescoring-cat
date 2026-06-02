## ============================================================
## SENSITIVITY CHECK (Ben comment #1): Wider difficulty spread,
## MULTI-POOL version.
##
## Re-runs the cross-population family with item difficulties
## drawn from a WIDER distribution, b ~ N(0, 1.5^2) instead of
## b ~ N(0, 0.75^2), across MULTIPLE independent item pools so
## the reported numbers are stable rather than a single-draw
## artifact.
##
## Why multi-pool: a single wide-spread pool produces a jumpy
## per-cell bias maximum (the metric that matters most here is a
## ratio, bias^2/MSE, taken as a max across 20 rho cells). With
## one pool at 50 reps that max is noisy. Averaging across pools,
## and reporting the distribution of the per-pool max, gives a
## number we can defend. This mirrors the logic of the paper's
## existing pool-robustness appendix.
##
## SELF-CONTAINED. Does NOT source any other script. Does NOT
## read or overwrite output/data/res.rds or any main-sim output.
##
## Settings match Study 1 EXACTLY except:
##   - b_sd:   0.75 -> 1.5         (the assumption under test)
##   - n_pools: 1   -> 10          (for stable estimates)
## Within each pool: 100 items, 100 examinees, 50 reps,
## 20 rho levels, 20-item CAT, EAP N(0,1.5) prior.
##
## Runtime: ~10x the single-pool run, roughly 4-7 hours.
## Run it unattended (e.g. nohup / overnight).
##
## USAGE: Rscript run_widespread_sensitivity_multipool.R
## ============================================================

library(tidyverse)

# ── Settings (match main sim EXCEPT b_sd and n_pools) ─────────
n_items     <- 100
n_examinees <- 100
n_reps      <- 50
cat_length  <- 20
prior_sd    <- 1.5
grid        <- seq(-6, 6, by = 0.05)

b_sd_wide   <- 1.5     # <<< assumption under test (main sim: 0.75)
n_pools     <- 10      # <<< independent pools for stability

rho_levels  <- seq(0.01, 0.99, length.out = 20)

out_dir <- "output/data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Generation (b_sd parameterised) ───────────────────────────
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

# ── Variance-preserving cross-population misalignment ─────────
# mu and sd_b are read from the pool, so a wider pool flows
# through cleanly; nothing re-centers or re-standardizes.
add_correlated <- function(pool, rho) {
  pn <- pool; mu <- mean(pool$b); sd_b <- sd(pool$b)
  residual_sd <- sd_b * sqrt(max(1 - rho^2, 0))
  pn$b <- mu + rho * (pool$b - mu) + rnorm(nrow(pool), 0, residual_sd)
  pn
}

# ── Evaluate one cell (per-rep metrics) ───────────────────────
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

# ── Per-pool summary: the two numbers, computed within a pool ──
summarise_pool <- function(res_pool) {
  res_cell <- res_pool %>%
    group_by(rho, condition) %>%
    summarise(
      cell_rmse     = sqrt(mean(mse)),
      cell_pct_bias = if_else(mean(mse) > 1e-10, mean(bias)^2 / mean(mse) * 100, 0),
      .groups = "drop"
    )
  
  scoring <- res_cell %>%
    select(rho, condition, cell_rmse) %>%
    pivot_wider(names_from = condition, values_from = cell_rmse) %>%
    mutate(
      total_error   = `CAT-O` - `CAT-B`,
      scoring_error = `CAT-O` - `CAT-R`,
      pct_scoring   = if_else(total_error > 1e-6, scoring_error / total_error * 100, NA_real_)
    )
  
  catr <- res_cell %>% filter(condition == "CAT-R")
  cato <- res_cell %>% filter(condition == "CAT-O")
  
  tibble(
    mean_pct_scoring = mean(scoring$pct_scoring, na.rm = TRUE),
    catr_bias_mean   = mean(catr$cell_pct_bias),
    catr_bias_max    = max(catr$cell_pct_bias),
    cato_bias_max    = max(cato$cell_pct_bias)
  )
}

# ══════════════════════════════════════════════════════════════
# RUN: 10 independent wide-spread pools
# ══════════════════════════════════════════════════════════════
set.seed(2026)  # reproducible: seeds the whole sequence of pools

pool_summaries <- vector("list", n_pools)
t0 <- Sys.time()

for (pp in seq_len(n_pools)) {
  cat(sprintf("\n=== POOL %d / %d  [%s] ===\n", pp, n_pools, format(Sys.time() - t0)))
  pool <- generate_pool(n_items, b_sd_wide)
  cat(sprintf("  SD(b) = %.3f | range [%.2f, %.2f]\n",
              sd(pool$b), min(pool$b), max(pool$b)))
  
  pool_res <- vector("list", length(rho_levels))
  for (r in seq_along(rho_levels)) {
    rho <- rho_levels[r]
    mis  <- add_correlated(pool, rho)
    cell <- evaluate_cell(pool, mis, n_examinees, n_reps, cat_length, grid, prior_sd)
    pool_res[[r]] <- compute_metrics_by_rep(cell) %>% mutate(rho = rho)
  }
  s <- summarise_pool(bind_rows(pool_res)) %>% mutate(pool = pp)
  pool_summaries[[pp]] <- s
  cat(sprintf("  pool scoring%%=%.1f | CAT-R bias mean=%.4f max=%.4f | CAT-O bias max=%.2f\n",
              s$mean_pct_scoring, s$catr_bias_mean, s$catr_bias_max, s$cato_bias_max))
}

pool_tab <- bind_rows(pool_summaries)

# ══════════════════════════════════════════════════════════════
# AGGREGATE ACROSS POOLS  (these are the numbers for the paper)
# ══════════════════════════════════════════════════════════════
cat("\n\n========== MULTI-POOL WIDE-SPREAD RESULTS",
    sprintf("(b_sd = %.2f, %d pools) ==========\n", b_sd_wide, n_pools))

cat("\nPer-pool values:\n")
print(pool_tab %>% select(pool, mean_pct_scoring, catr_bias_mean,
                          catr_bias_max, cato_bias_max))

agg <- pool_tab %>%
  summarise(
    scoring_mean   = mean(mean_pct_scoring),
    scoring_min    = min(mean_pct_scoring),
    scoring_max    = max(mean_pct_scoring),
    catr_bias_mean = mean(catr_bias_mean),
    # mean of the per-pool maxima, and the worst single cell across all pools
    catr_max_meanofmax = mean(catr_bias_max),
    catr_max_worst     = max(catr_bias_max),
    cato_bias_max  = max(cato_bias_max)
  )

cat("\n--- Aggregated across pools (report these) ---\n")
cat(sprintf("  Scoring error share: mean %.1f%%  (range across pools %.1f%% to %.1f%%)\n",
            agg$scoring_mean, agg$scoring_min, agg$scoring_max))
cat(sprintf("  CAT-R bias: mean %.4f%% of MSE\n", agg$catr_bias_mean))
cat(sprintf("  CAT-R bias per-pool max: mean of maxima %.4f%%, worst single cell %.4f%%\n",
            agg$catr_max_meanofmax, agg$catr_max_worst))
cat(sprintf("  CAT-O bias max (contrast): %.2f%%\n", agg$cato_bias_max))

cat("\n--- Interpretation guide ---\n")
cat("  Report scoring_mean and its across-pool range as the headline.\n")
cat("  For bias, the defensible statement is the across-pool MEAN and the\n")
cat("  mean-of-per-pool-maxima. The single worst cell is the noisiest number;\n")
cat("  do not lead with it.\n")

write_csv(pool_tab, file.path(out_dir, "sensitivity_widespread_multipool.csv"))
cat("\nSaved: output/data/sensitivity_widespread_multipool.csv\n")
cat("\nDone!\n")