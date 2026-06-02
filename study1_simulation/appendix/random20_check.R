## ============================================================
## Standalone check: Random-20 vs CAT-B vs Full-100
## Answers Ben's comment: how does scoring 20 RANDOM items
## compare to adaptive-20 (CAT-B) and the full 100-item ceiling?
##
## Self-contained: defines its own functions, does NOT source
## run_full_simulation_v2.R, does NOT touch output/data/res.rds.
##
## Random-20 RMSE is test-length cost and is independent of rho,
## so a single representative cell suffices.
##
## Runtime: ~5 min
## USAGE: Rscript run_random20_check.R
## ============================================================

library(tidyverse)

# ── Settings (match main sim exactly) ─────────────────────────
n_items     <- 100
n_examinees <- 100
n_reps      <- 50
cat_length  <- 20
prior_sd    <- 1.5
grid        <- seq(-6, 6, by = 0.05)
rho_test    <- 0.5    # any value; Random-20 does not depend on rho

# ── Generation ────────────────────────────────────────────────
generate_pool <- function(n_items) {
  tibble(
    item = 1:n_items,
    b = rnorm(n_items, 0, 0.75),
    a = pmax(0.5, exp(rnorm(n_items, 0.1, 0.3)))
  )
}
generate_examinees <- function(n) rnorm(n, 0, 1)

# ── IRT ───────────────────────────────────────────────────────
p_2pl <- function(theta, a, b) 1 / (1 + exp(-a * (theta - b)))
info_2pl <- function(theta, a, b) {
  p <- p_2pl(theta, a, b); a^2 * p * (1 - p)
}
eap_score <- function(responses, a, b, grid, prior_sd) {
  ll <- sapply(grid, function(th) {
    p <- p_2pl(th, a, b); p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
    sum(responses * log(p) + (1 - responses) * log(1 - p))
  })
  lp <- ll + dnorm(grid, 0, prior_sd, log = TRUE)
  w  <- exp(lp - max(lp)); w <- w / sum(w)
  sum(grid * w)
}

# ── CAT + fixed-set scoring ───────────────────────────────────
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
  list(theta = theta_final, items = items_given)
}
run_fixed <- function(theta_true, pool, items_given, score_bank, grid, prior_sd) {
  responses <- rbinom(length(items_given), 1,
                      p_2pl(theta_true, pool$a[items_given], pool$b[items_given]))
  theta <- eap_score(responses, score_bank$a[items_given],
                     score_bank$b[items_given], grid, prior_sd)
  list(theta = theta)
}
add_correlated <- function(pool, rho) {
  pn <- pool; mu <- mean(pool$b); sd_b <- sd(pool$b)
  residual_sd <- sd_b * sqrt(max(1 - rho^2, 0))
  pn$b <- mu + rho * (pool$b - mu) + rnorm(nrow(pool), 0, residual_sd)
  pn
}

# ── Run one cell, all 6 conditions ────────────────────────────
set.seed(2026)
pool <- generate_pool(n_items)
mis  <- add_correlated(pool, rho_test)
true_bank  <- list(a = pool$a, b = pool$b)
noisy_bank <- list(a = mis$a,  b = mis$b)
fixed_items <- order(pool$a, decreasing = TRUE)[1:cat_length]

rows <- list(); idx <- 0
for (rep in seq_len(n_reps)) {
  thetas <- generate_examinees(n_examinees)
  for (i in seq_along(thetas)) {
    th <- thetas[i]
    cat_o <- run_cat(th, pool, noisy_bank, noisy_bank, cat_length, grid, prior_sd)
    cat_b <- run_cat(th, pool, true_bank, true_bank, cat_length, grid, prior_sd)
    fix   <- run_fixed(th, pool, fixed_items, true_bank, grid, prior_sd)
    full  <- run_fixed(th, pool, seq_len(n_items), true_bank, grid, prior_sd)
    rand_items <- sample(seq_len(n_items), cat_length)
    randd <- run_fixed(th, pool, rand_items, true_bank, grid, prior_sd)
    b <- tibble(theta_true = th)
    idx <- idx+1; rows[[idx]] <- b %>% mutate(condition="CAT-O",     theta_hat=cat_o$theta)
    idx <- idx+1; rows[[idx]] <- b %>% mutate(condition="CAT-B",     theta_hat=cat_b$theta)
    idx <- idx+1; rows[[idx]] <- b %>% mutate(condition="Fixed-20",  theta_hat=fix$theta)
    idx <- idx+1; rows[[idx]] <- b %>% mutate(condition="Full-100",  theta_hat=full$theta)
    idx <- idx+1; rows[[idx]] <- b %>% mutate(condition="Random-20", theta_hat=randd$theta)
  }
}
cell <- bind_rows(rows)

cat("\n========== RMSE by condition (rho =", rho_test, ") ==========\n")
cell %>%
  group_by(condition) %>%
  summarise(rmse = sqrt(mean((theta_hat - theta_true)^2)), .groups="drop") %>%
  arrange(rmse) %>%
  print()
cat("\nKey comparison for Ben: Random-20 vs CAT-B (adaptive 20) vs Full-100.\n")
cat("If CAT-B < Random-20, adaptive selection helps beyond just test length.\n")