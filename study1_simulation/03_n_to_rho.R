## ============================================================
## Secondary Simulation: N → ρ Mapping
## Maps respondent sample size to item parameter recovery (ρ)
## Matches main simulation settings from run_main_simulation.R
##
## Output: output/data/n_to_rho.rds, figures/paper/fig_n_to_rho.png
## Runtime: ~10 min
## ============================================================

library(tidyverse)
library(mirt)

library(here)   # paths resolve relative to repo root (.here marker)
set.seed(2026)

# ── Settings (match main simulation) ──────────────────────────

n_items   <- 100
b_sd      <- 0.75       # same as main sim: b ~ N(0, 0.75^2)
n_reps    <- 200        # more reps since this is cheap

N_vec <- c(5, 8, 10, 15, 20, 25, 50, 75, 100, 150, 200, 300, 500)
# ── Run ───────────────────────────────────────────────────────

cat("Running N-to-rho simulation...\n")
cat("Items:", n_items, "| b_sd:", b_sd, "| Reps:", n_reps, "\n")
cat("N values:", paste(N_vec, collapse = ", "), "\n\n")

results <- expand_grid(N = N_vec, rep = 1:n_reps)

# Use a FIXED item pool across all reps (cleaner for showing N effect)
b_true <- rnorm(n_items, 0, b_sd)

results$rho <- NA_real_

pb <- txtProgressBar(min = 0, max = nrow(results), style = 3)
for (i in seq_len(nrow(results))) {
  n <- results$N[i]
  
  # Generate responses (Rasch / 1PL)
  theta <- rnorm(n, 0, 1)
  prob  <- plogis(outer(theta, b_true, "-"))
  resp  <- matrix(rbinom(n * n_items, 1, prob), n, n_items)
  colnames(resp) <- paste0("i", 1:n_items)
  
  # Estimate difficulties via marginal MLE (Rasch)
  tryCatch({
    mod <- mirt(as.data.frame(resp), 1, itemtype = "Rasch",
                verbose = FALSE, technical = list(NCYCLES = 500))
    # Extract difficulty: in mirt, d parameter = -b for Rasch
    coefs <- coef(mod, simplify = TRUE)$items
    b_hat <- -coefs[, "d"]
    
    results$rho[i] <- cor(b_true, b_hat)
  }, error = function(e) {
    results$rho[i] <<- NA
  })
  
  setTxtProgressBar(pb, i)
}
close(pb)

# ── Summarise ─────────────────────────────────────────────────

n_to_rho <- results %>%
  filter(!is.na(rho)) %>%
  group_by(N) %>%
  summarise(
    rho_mean = mean(rho),
    rho_sd   = sd(rho),
    rho_se   = sd(rho) / sqrt(n()),
    rho_q25  = quantile(rho, 0.25),
    rho_q75  = quantile(rho, 0.75),
    n_valid  = n(),
    .groups  = "drop"
  )

cat("\n\n========== N → ρ MAPPING ==========\n")
print(n_to_rho)

dir.create(here("output", "data"), showWarnings = FALSE, recursive = TRUE)
saveRDS(results,  here("output", "data", "n_to_rho_raw.rds"))
saveRDS(n_to_rho, here("output", "data", "n_to_rho.rds"))
cat("\nSaved output/data/n_to_rho.rds and n_to_rho_raw.rds\n")

results  <- readRDS(here("output", "data", "n_to_rho_raw.rds"))
n_to_rho <- readRDS(here("output", "data", "n_to_rho.rds"))

# ── Plot ──────────────────────────────────────────────────────

p <- ggplot(n_to_rho, aes(x = N, y = rho_mean)) +
  geom_ribbon(aes(ymin = rho_q25, ymax = rho_q75),
              fill = "#E69F00", alpha = 0.2) +
  geom_line(color = "#E69F00", linewidth = 1.2) +
  geom_point(color = "#E69F00", size = 2.5) +
  geom_errorbar(aes(ymin = rho_mean - rho_se, ymax = rho_mean + rho_se),
                width = 15, color = "#E69F00", alpha = 0.6) +
  # Reference lines for LEVANTE empirical values
  geom_hline(yintercept = 0.34, linetype = "dashed", color = "#D55E00", alpha = 0.6) +
  annotate("text", x = 900, y = 0.36, label = "Matrix reasoning (ρ = .34)",
           size = 3, color = "#D55E00", family = "serif") +
  geom_hline(yintercept = 0.83, linetype = "dashed", color = "#0072B2", alpha = 0.6) +
  annotate("text", x = 900, y = 0.85, label = "Same-different (ρ = .83)",
           size = 3, color = "#0072B2", family = "serif") +
  scale_x_continuous(breaks = N_vec) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
  labs(
    x = "Respondent sample size (N)",
    y = expression(Item~parameter~recovery~(rho))
  ) +
  theme_minimal(base_size = 13) +
  theme(
    text            = element_text(family = "serif"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#eeeeee", linewidth = 0.3),
    axis.title      = element_text(size = 12),
    axis.text       = element_text(size = 9, color = "#444444")
  )

dir.create(here("figures", "paper"), showWarnings = FALSE, recursive = TRUE)
ggsave(here("figures", "paper", "fig_n_to_rho.png"), p,
       width = 8, height = 4.5, dpi = 300, bg = "white")
cat("Saved figures/paper/fig_n_to_rho.png\n")

# ── Print interpretive summary ────────────────────────────────

cat("\n========== INTERPRETATION ==========\n")
cat("At LEVANTE's empirical ρ values:\n")
for (target_rho in c(0.34, 0.41, 0.70, 0.83)) {
  # Interpolate N from the mapping
  if (target_rho < min(n_to_rho$rho_mean)) {
    cat(sprintf("  ρ = %.2f → below N = %d (ρ_mean = %.2f)\n",
                target_rho, min(N_vec), min(n_to_rho$rho_mean)))
  } else {
    fit <- approx(n_to_rho$rho_mean, n_to_rho$N, xout = target_rho)
    cat(sprintf("  ρ = %.2f → equivalent to N ≈ %d respondents\n",
                target_rho, round(fit$y)))
  }
}

cat("\nDone!\n")

