## ============================================================
## Full CAT Rescoring Simulation
##
## Three error families × 5 conditions + stepwise trajectory
## Outputs: res.rds, res_stepwise.rds
##
## Run from: ~/Documents/levante-runs/levante-runs/rescoring_simulation
## Runtime: ~3-4 hours. Quick test: set n_reps=5, n_examinees=20
## ============================================================

library(tidyverse)

set.seed(2026)

# ── 0. Settings (match paper exactly) ─────────────────────────

n_items      <- 100       # 2PL pool
n_examinees  <- 100       # per cell
n_reps       <- 50        # replications per cell
cat_length   <- 20        # fixed CAT length
prior_sd     <- 1.5       # EAP prior N(0, 1.5)
grid         <- seq(-6, 6, by = 0.05)

out_dir <- "output/data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ── 1. Item & examinee generation ─────────────────────────────
# Paper: b ~ N(0, 0.75^2), a ~ lognormal(0.1, 0.3) floor 0.5

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


# ── 3. CAT engine (with optional trajectory) ──────────────────

run_cat <- function(theta_true, pool, select_bank, score_bank,
                    cat_length, grid, prior_sd,
                    return_trajectory = FALSE) {
  n_pool    <- nrow(pool)
  used      <- logical(n_pool)
  responses <- numeric(0)
  items_given <- integer(0)
  theta_hat <- 0
  
  # For trajectory: score at each step with BOTH banks
  traj_score <- traj_select <- numeric(cat_length)
  
  for (k in seq_len(cat_length)) {
    # Select using select_bank
    info_vals <- info_2pl(theta_hat, select_bank$a, select_bank$b)
    info_vals[used] <- -Inf
    chosen <- which.max(info_vals)
    used[chosen] <- TRUE
    items_given <- c(items_given, chosen)
    
    # Generate response using TRUE parameters
    p <- p_2pl(theta_true, pool$a[chosen], pool$b[chosen])
    resp <- rbinom(1, 1, p)
    responses <- c(responses, resp)
    
    # Update theta for next selection (using select_bank)
    theta_hat <- eap_score(responses, select_bank$a[items_given],
                           select_bank$b[items_given], grid, prior_sd)
    
    if (return_trajectory) {
      traj_select[k] <- theta_hat
      traj_score[k]  <- eap_score(responses, score_bank$a[items_given],
                                  score_bank$b[items_given], grid, prior_sd)
    }
  }
  
  # Final scoring with score_bank
  theta_final <- eap_score(responses, score_bank$a[items_given],
                           score_bank$b[items_given], grid, prior_sd)
  
  out <- list(theta = theta_final, items = items_given, responses = responses)
  if (return_trajectory) {
    out$traj_select <- traj_select
    out$traj_score  <- traj_score
  }
  out
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

# Random noise: b' = b + N(0, sigma^2)
add_noise <- function(pool, sigma) {
  pool_noisy <- pool
  if (sigma == 0) return(pool_noisy)  # no noise
  pool_noisy$b <- pool$b + rnorm(nrow(pool), 0, sigma)
  pool_noisy
}

# Structured shift + noise: b' = b + delta + N(0, sigma^2)
add_shift_noise <- function(pool, delta, sigma) {
  pool_noisy <- pool
  pool_noisy$b <- pool$b + delta + rnorm(nrow(pool), 0, sigma)
  pool_noisy
}

# Cross-population: b' = mu + rho*(b - mu) + epsilon
# epsilon ~ N(0, sd_b * sqrt(1 - rho^2)) so cor(b, b') ≈ rho
add_correlated <- function(pool, rho) {
  pool_noisy <- pool
  mu <- mean(pool$b)
  sd_b <- sd(pool$b)
  residual_sd <- sd_b * sqrt(max(1 - rho^2, 0))
  pool_noisy$b <- mu + rho * (pool$b - mu) + rnorm(nrow(pool), 0, residual_sd)
  pool_noisy
}


# ── 5. Evaluate one cell ──────────────────────────────────────

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
      
      # CAT-O: select noisy, score noisy
      cat_o <- run_cat(th, pool, noisy_bank, noisy_bank,
                       cat_length, grid, prior_sd)
      
      # CAT-R: same items as CAT-O, rescore with true
      theta_r <- eap_score(cat_o$responses, true_bank$a[cat_o$items],
                           true_bank$b[cat_o$items], grid, prior_sd)
      
      # CAT-B: select true, score true
      cat_b <- run_cat(th, pool, true_bank, true_bank,
                       cat_length, grid, prior_sd)
      
      # Fixed-20
      fix <- run_fixed(th, pool, fixed_items, true_bank, grid, prior_sd)
      
      # Full-100
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


# ── 6. Compute metrics ────────────────────────────────────────

compute_metrics <- function(cell_results) {
  cell_results %>%
    group_by(condition) %>%
    summarise(
      rmse     = sqrt(mean((theta_hat - theta_true)^2)),
      bias     = mean(theta_hat - theta_true),
      cor_r    = cor(theta_hat, theta_true),
      mse      = mean((theta_hat - theta_true)^2),
      bias_sq  = mean(theta_hat - theta_true)^2,
      variance = mse - bias_sq,
      pct_bias = if_else(mse > 1e-10, bias_sq / mse * 100, 0),
      n        = n(),
      .groups  = "drop"
    )
}


# ══════════════════════════════════════════════════════════════
# 7. RUN ALL THREE FAMILIES
# ══════════════════════════════════════════════════════════════

all_results <- list()
t0 <- Sys.time()

# --- Family 1: Cross-population misalignment (rho) ---
cat("=== Family 1: Cross-population misalignment ===\n")
rho_levels <- seq(0.01, 0.99, length.out = 20)

for (r in seq_along(rho_levels)) {
  rho <- rho_levels[r]
  cat(sprintf("  rho = %.3f (%d/%d) [%s]\n", rho, r, length(rho_levels),
              format(Sys.time() - t0)))
  
  pool <- generate_pool(n_items)
  misaligned <- add_correlated(pool, rho)
  
  cell <- evaluate_cell(pool, misaligned, n_examinees, n_reps,
                        cat_length, grid, prior_sd)
  metrics <- compute_metrics(cell)
  
  all_results <- c(all_results, list(
    metrics %>% mutate(mode = "correlated", rho = rho,
                       sigma = NA_real_, delta = NA_real_)
  ))
}


# --- Family 2: Random noise (sigma) ---
# Paper: sigma in [0, 1], 20 levels. Include sigma=0 as baseline.
cat("\n=== Family 2: Random noise ===\n")
sigma_levels <- c(0, seq(0.05, 1.0, length.out = 19))  # 20 levels, starting at 0

for (s in seq_along(sigma_levels)) {
  sig <- sigma_levels[s]
  cat(sprintf("  sigma = %.3f (%d/%d) [%s]\n", sig, s, length(sigma_levels),
              format(Sys.time() - t0)))
  
  pool <- generate_pool(n_items)
  misaligned <- add_noise(pool, sig)
  
  cell <- evaluate_cell(pool, misaligned, n_examinees, n_reps,
                        cat_length, grid, prior_sd)
  metrics <- compute_metrics(cell)
  
  all_results <- c(all_results, list(
    metrics %>% mutate(mode = "random_noise", rho = NA_real_,
                       sigma = sig, delta = NA_real_)
  ))
}


# --- Family 3: Structured shift + noise ---
# Paper: delta in {-0.3, +0.3, +1.0}, sigma varies
cat("\n=== Family 3: Structured shifts ===\n")
delta_levels <- c(-0.3, 0.3, 1.0)
sigma_shift_levels <- seq(0.05, 1.0, length.out = 10)

for (d in delta_levels) {
  for (s in seq_along(sigma_shift_levels)) {
    sig <- sigma_shift_levels[s]
    cat(sprintf("  delta = %.1f, sigma = %.3f (%d/%d) [%s]\n",
                d, sig, s, length(sigma_shift_levels),
                format(Sys.time() - t0)))
    
    pool <- generate_pool(n_items)
    misaligned <- add_shift_noise(pool, d, sig)
    
    cell <- evaluate_cell(pool, misaligned, n_examinees, n_reps,
                          cat_length, grid, prior_sd)
    metrics <- compute_metrics(cell)
    
    all_results <- c(all_results, list(
      metrics %>% mutate(mode = "structured_shift", rho = NA_real_,
                         sigma = sig, delta = d)
    ))
  }
}


# ══════════════════════════════════════════════════════════════
# 8. STEPWISE TRAJECTORY (at rho = 0.75)
# ══════════════════════════════════════════════════════════════

cat("\n=== Stepwise trajectory (rho = 0.75, n = 300) ===\n")

n_step_examinees <- 300
n_step_reps      <- 50
rho_step         <- 0.75

pool_step      <- generate_pool(n_items)
misaligned_step <- add_correlated(pool_step, rho_step)

true_bank_s  <- list(a = pool_step$a, b = pool_step$b)
noisy_bank_s <- list(a = misaligned_step$a, b = misaligned_step$b)
fixed_items_s <- order(pool_step$a, decreasing = TRUE)  # all items ranked

step_results <- list()

for (rep in seq_len(n_step_reps)) {
  if (rep %% 10 == 0) cat(sprintf("  rep %d/%d [%s]\n", rep, n_step_reps,
                                   format(Sys.time() - t0)))
  thetas <- generate_examinees(n_step_examinees)
  
  for (i in seq_along(thetas)) {
    th <- thetas[i]
    
    # CAT-O with trajectory (score with noisy)
    cat_o <- run_cat(th, pool_step, noisy_bank_s, noisy_bank_s,
                     cat_length, grid, prior_sd, return_trajectory = TRUE)
    
    # CAT-R trajectory: rescore CAT-O items at each step with true bank
    traj_r <- numeric(cat_length)
    for (k in seq_len(cat_length)) {
      items_k <- cat_o$items[1:k]
      resp_k  <- cat_o$responses[1:k]
      traj_r[k] <- eap_score(resp_k, true_bank_s$a[items_k],
                             true_bank_s$b[items_k], grid, prior_sd)
    }
    
    # CAT-B with trajectory
    cat_b <- run_cat(th, pool_step, true_bank_s, true_bank_s,
                     cat_length, grid, prior_sd, return_trajectory = TRUE)
    
    # Fixed-form trajectory (first k of ranked items)
    traj_fix <- numeric(cat_length)
    all_fix_resp <- rbinom(cat_length, 1,
                           p_2pl(th, pool_step$a[fixed_items_s[1:cat_length]],
                                 pool_step$b[fixed_items_s[1:cat_length]]))
    for (k in seq_len(cat_length)) {
      items_k <- fixed_items_s[1:k]
      traj_fix[k] <- eap_score(all_fix_resp[1:k], true_bank_s$a[items_k],
                               true_bank_s$b[items_k], grid, prior_sd)
    }
    
    base <- tibble(rep = rep, examinee = i, theta_true = th)
    
    for (k in seq_len(cat_length)) {
      step_results <- c(step_results, list(
        base %>% mutate(step = k,
                        CAT_O = cat_o$traj_select[k],  # scored with noisy
                        CAT_R = traj_r[k],
                        CAT_B = cat_b$traj_score[k],
                        Fixed = traj_fix[k])
      ))
    }
  }
}

stepwise <- bind_rows(step_results) %>%
  pivot_longer(cols = c(CAT_O, CAT_R, CAT_B, Fixed),
               names_to = "condition", values_to = "theta_hat") %>%
  mutate(condition = recode(condition,
                            "CAT_O" = "CAT-O", "CAT_R" = "CAT-R",
                            "CAT_B" = "CAT-B", "Fixed" = "Fixed-20")) %>%
  group_by(condition, step) %>%
  summarise(
    rmse = sqrt(mean((theta_hat - theta_true)^2)),
    bias = mean(theta_hat - theta_true),
    cor_r = cor(theta_hat, theta_true),
    .groups = "drop"
  )


# ══════════════════════════════════════════════════════════════
# 9. SAVE
# ══════════════════════════════════════════════════════════════

res <- bind_rows(all_results)

cat("\n\n========== SUMMARY ==========\n")
cat("Total time:", format(Sys.time() - t0), "\n")
cat("res.rds:", nrow(res), "rows\n")
count(res, mode) %>% print()
cat("\nstepwise:", nrow(stepwise), "rows\n")

saveRDS(res, file.path(out_dir, "res.rds"))
saveRDS(stepwise, file.path(out_dir, "res_stepwise.rds"))

cat("\nSaved:\n")
cat("  ", file.path(out_dir, "res.rds"), "\n")
cat("  ", file.path(out_dir, "res_stepwise.rds"), "\n")


# ══════════════════════════════════════════════════════════════
# 10. QUICK CHECKS
# ══════════════════════════════════════════════════════════════

cat("\n-- Correlated: CAT-O vs CAT-R at low rho --\n")
res %>%
  filter(mode == "correlated", condition %in% c("CAT-O", "CAT-R")) %>%
  filter(rho < 0.15) %>%
  select(condition, rho, rmse, bias, pct_bias) %>%
  print()

cat("\n-- Random noise: sigma = 0 baseline --\n")
res %>%
  filter(mode == "random_noise", sigma == 0) %>%
  select(condition, sigma, rmse, bias) %>%
  print()

cat("\n-- Bias < 0.1% of MSE across all families? --\n")
res %>%
  group_by(mode) %>%
  summarise(max_pct_bias = max(pct_bias), .groups = "drop") %>%
  print()

cat("\n-- Scoring error > selection error? --\n")
res %>%
  select(mode, rho, sigma, delta, condition, rmse) %>%
  pivot_wider(names_from = condition, values_from = rmse) %>%
  mutate(
    scoring_err   = `CAT-O` - `CAT-R`,
    selection_err = `CAT-R` - `CAT-B`,
    pct_scoring   = scoring_err / (`CAT-O` - `CAT-B`) * 100
  ) %>%
  group_by(mode) %>%
  summarise(mean_pct_scoring = mean(pct_scoring, na.rm = TRUE),
            .groups = "drop") %>%
  print()

cat("\n-- Stepwise: RMSE at steps 5, 10, 15, 20 --\n")
stepwise %>%
  filter(step %in% c(5, 10, 15, 20)) %>%
  arrange(step, condition) %>%
  print(n = 20)

cat("\nDone!\n")
