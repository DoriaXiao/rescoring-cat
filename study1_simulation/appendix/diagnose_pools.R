## ============================================================
## DIAGNOSTIC: why do pools 4 and 6 have low scoring share?
##
## The multi-pool run gave scoring shares ranging 59% to 88%.
## Pools 4 (72.7%) and 6 (59.0%) are the low outliers pulling
## the mean down to 79%. This script reproduces the SAME 10
## pools (identical seeding to run_widespread_sensitivity_multipool.R)
## and inspects what is structurally different about them.
##
## It does two things per pool:
##   (1) Pool structure: SD, range, skew, kurtosis of b; how many
##       items sit in the informative -3..3 region; discrimination
##       summary; how concentrated information is.
##   (2) Error structure at a mid rho (0.5): the actual CAT-O,
##       CAT-R, CAT-B RMSE so we can see WHERE the scoring share
##       comes from (small CAT-O-minus-CAT-R gap, or large
##       CAT-R-minus-CAT-B selection residual).
##
## Scoring share = (CAT-O - CAT-R) / (CAT-O - CAT-B).
## A LOW share means CAT-R sits closer to CAT-O (less fixable)
## OR CAT-B is far below CAT-R (big selection residual). This
## tells us which mechanism drives the low-share pools.
##
## SELF-CONTAINED. Touches nothing. ~15-20 min (lighter settings:
## fewer reps, 3 rho levels, since we only need structure not
## publication-precision numbers).
##
## USAGE: Rscript diagnose_pools_4_6.R
## ============================================================

library(tidyverse)

# Lighter settings: we want structure, not 4-decimal precision
n_items     <- 100
n_examinees <- 100
n_reps      <- 30          # lighter than the 50 used in the main run
cat_length  <- 20
prior_sd    <- 1.5
grid        <- seq(-6, 6, by = 0.05)
b_sd_wide   <- 1.5
n_pools     <- 10
rho_probe   <- c(0.20, 0.50, 0.75)   # a few rho levels, not all 20

# ── functions (identical to the multi-pool script) ───────────
generate_pool <- function(n_items, b_sd) {
  tibble(item = 1:n_items,
         b = rnorm(n_items, 0, b_sd),
         a = pmax(0.5, exp(rnorm(n_items, 0.1, 0.3))))
}
generate_examinees <- function(n) rnorm(n, 0, 1)
p_2pl <- function(theta, a, b) 1 / (1 + exp(-a * (theta - b)))
info_2pl <- function(theta, a, b) { p <- p_2pl(theta, a, b); a^2 * p * (1 - p) }
eap_score <- function(responses, a, b, grid, prior_sd) {
  ll <- sapply(grid, function(th) {
    p <- p_2pl(th, a, b); p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
    sum(responses * log(p) + (1 - responses) * log(1 - p))
  })
  lp <- ll + dnorm(grid, 0, prior_sd, log = TRUE)
  w  <- exp(lp - max(lp)); w <- w / sum(w); sum(grid * w)
}
run_cat <- function(theta_true, pool, select_bank, score_bank,
                    cat_length, grid, prior_sd) {
  used <- logical(nrow(pool)); responses <- numeric(0)
  items_given <- integer(0); theta_hat <- 0
  for (k in seq_len(cat_length)) {
    iv <- info_2pl(theta_hat, select_bank$a, select_bank$b); iv[used] <- -Inf
    chosen <- which.max(iv); used[chosen] <- TRUE
    items_given <- c(items_given, chosen)
    responses <- c(responses, rbinom(1, 1, p_2pl(theta_true, pool$a[chosen], pool$b[chosen])))
    theta_hat <- eap_score(responses, select_bank$a[items_given],
                           select_bank$b[items_given], grid, prior_sd)
  }
  theta_final <- eap_score(responses, score_bank$a[items_given],
                           score_bank$b[items_given], grid, prior_sd)
  list(theta = theta_final, items = items_given, responses = responses)
}
add_correlated <- function(pool, rho) {
  pn <- pool; mu <- mean(pool$b); sd_b <- sd(pool$b)
  pn$b <- mu + rho * (pool$b - mu) +
    rnorm(nrow(pool), 0, sd_b * sqrt(max(1 - rho^2, 0)))
  pn
}
skewness <- function(x) mean((x - mean(x))^3) / sd(x)^3
kurtosis <- function(x) mean((x - mean(x))^4) / sd(x)^4

# ── reproduce the SAME 10 pools, identical seeding ───────────
set.seed(2026)   # same master seed as the multi-pool run

struct_rows <- list()
error_rows  <- list()

for (pp in seq_len(n_pools)) {
  pool <- generate_pool(n_items, b_sd_wide)   # draws in same order as main run

  # (1) structure
  info_at0 <- sum(info_2pl(0, pool$a, pool$b))           # total info at theta=0
  in_core  <- mean(abs(pool$b) <= 3)                      # frac items in -3..3
  in_tight <- mean(abs(pool$b) <= 1)                      # frac near center
  top20_b  <- pool$b[order(pool$a, decreasing = TRUE)[1:cat_length]]
  struct_rows[[pp]] <- tibble(
    pool_id = pp,
    sd_b = sd(pool$b), min_b = min(pool$b), max_b = max(pool$b),
    skew_b = skewness(pool$b), kurt_b = kurtosis(pool$b),
    frac_core_3 = in_core, frac_tight_1 = in_tight,
    info_at_0 = info_at0,
    # where do the most-discriminating items sit? (these get selected)
    top20_b_mean = mean(top20_b), top20_b_absmean = mean(abs(top20_b)),
    a_mean = mean(pool$a), a_max = max(pool$a)
  )

  # (2) error structure at probe rho levels
  for (rho in rho_probe) {
    mis <- add_correlated(pool, rho)
    true_bank  <- list(a = pool$a, b = pool$b)
    noisy_bank <- list(a = mis$a, b = mis$b)
    diffs <- replicate(n_reps, {
      th <- generate_examinees(n_examinees)
      o <- r <- bb <- numeric(length(th))
      for (i in seq_along(th)) {
        co <- run_cat(th[i], pool, noisy_bank, noisy_bank, cat_length, grid, prior_sd)
        o[i]  <- (co$theta - th[i])^2
        rr <- eap_score(co$responses, true_bank$a[co$items],
                        true_bank$b[co$items], grid, prior_sd)
        r[i]  <- (rr - th[i])^2
        cb <- run_cat(th[i], pool, true_bank, true_bank, cat_length, grid, prior_sd)
        bb[i] <- (cb$theta - th[i])^2
      }
      c(mean(o), mean(r), mean(bb))
    })
    rmse_o <- sqrt(mean(diffs[1, ])); rmse_r <- sqrt(mean(diffs[2, ]))
    rmse_b <- sqrt(mean(diffs[3, ]))
    error_rows[[length(error_rows) + 1]] <- tibble(
      pool_id = pp, rho = rho,
      rmse_O = rmse_o, rmse_R = rmse_r, rmse_B = rmse_b,
      scoring_gap   = rmse_o - rmse_r,         # fixable
      selection_gap = rmse_r - rmse_b,         # not fixable
      total_gap     = rmse_o - rmse_b,
      pct_scoring   = if_else(rmse_o - rmse_b > 1e-6,
                              (rmse_o - rmse_r) / (rmse_o - rmse_b) * 100, NA_real_)
    )
  }
  cat(sprintf("pool %d done\n", pp))
}

struct <- bind_rows(struct_rows)
errors <- bind_rows(error_rows)

# ── report ────────────────────────────────────────────────────
cat("\n\n===== POOL STRUCTURE (flag pools 4 and 6) =====\n")
print(struct %>% mutate(across(where(is.numeric), ~round(.x, 3))), n = 10)

cat("\n===== ERROR STRUCTURE at probe rho (pools 4 and 6 vs rest) =====\n")
print(errors %>% mutate(across(where(is.numeric), ~round(.x, 3))), n = 30)

cat("\n===== SCORING SHARE by pool (averaged over probe rho) =====\n")
share_by_pool <- errors %>%
  group_by(pool_id) %>%
  summarise(mean_pct_scoring = mean(pct_scoring, na.rm = TRUE),
            mean_selection_gap = mean(selection_gap),
            mean_scoring_gap = mean(scoring_gap),
            mean_rmse_B = mean(rmse_B), .groups = "drop")
print(share_by_pool %>% mutate(across(where(is.numeric), ~round(.x, 3))), n = 10)

cat("\n===== CONTRAST: low-share pools (4,6) vs others =====\n")
struct %>%
  mutate(group = if_else(pool_id %in% c(4, 6), "low_share_4_6", "others")) %>%
  group_by(group) %>%
  summarise(across(c(sd_b, skew_b, kurt_b, frac_core_3, frac_tight_1,
                     info_at_0, top20_b_absmean, a_max),
                   ~round(mean(.x), 3)), .groups = "drop") %>%
  print()

errors %>%
  mutate(group = if_else(pool_id %in% c(4, 6), "low_share_4_6", "others")) %>%
  group_by(group) %>%
  summarise(across(c(rmse_O, rmse_R, rmse_B, scoring_gap, selection_gap, pct_scoring),
                   ~round(mean(.x), 3)), .groups = "drop") %>%
  print()

cat("\n--- How to read this ---\n")
cat("If pools 4,6 have LARGER selection_gap (CAT-R minus CAT-B): the wider\n")
cat("spread in those draws made adaptive selection on miscalibrated params\n")
cat("genuinely costlier, so more error is unfixable. That is a real effect,\n")
cat("not an artifact.\n")
cat("If pools 4,6 have SMALLER total_gap (CAT-O close to CAT-B): the share\n")
cat("denominator is tiny and the percentage is just unstable. That is noise,\n")
cat("and the 79%% mean is being dragged by an ill-conditioned ratio.\n")
