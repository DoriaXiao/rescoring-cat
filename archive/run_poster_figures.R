## ============================================================
## Generate NCME Poster Figures 1 & 2
## Requires: res.rds (from run_main_simulation.R)
## Output:   poster_fig1_rho.png, poster_fig2_stepwise.png
## ============================================================

library(dplyr); library(ggplot2); library(purrr); library(tidyr)
library(ggrepel); library(withr); library(tibble)

# ── Source order matters ──────────────────────────────────────
# 1. simulation_source.R: core IRT functions (generate_items_pool, estimate_theta, etc.)
# 2. plot_traj.R: overrides add_calibration_error (err_mode) + plot_traj_simple (n_reps)
# 3. plotting_source.R: plot_metric_accuracy_v4 (does NOT call add_calibration_error)
source("R/simulation_source.R")
source("R/plot_traj.R")
source("R/plotting_source.R")

# ── Load simulation results ───────────────────────────────────
res <- readRDS("res.rds")
cat("Loaded res.rds:", nrow(res), "rows\n")
cat("Modes:", paste(unique(res$mode), collapse = ", "), "\n")

# ══════════════════════════════════════════════════════════════
# FIGURE 1: RMSE vs ρ (correlated misalignment)
# ══════════════════════════════════════════════════════════════
# This is the main conceptual figure for the poster.
# Shows: CAT-O degrades as ρ drops; CAT-R stays flat;
#        Fixed-20 is constant; gap = selection error.

fig1 <- plot_metric_accuracy_v4(
  res,
  metric  = "rmse",
  family  = "correlated",
  xvar    = "rho",
  smooth  = TRUE,
  smooth_method = "loess",
  smooth_span   = 0.5,
  error_bar = TRUE,
  show_labels = TRUE
)

# Override title/subtitle for poster
fig1 <- fig1 +
  labs(
    title    = "Performance vs Alignment ρ",
    subtitle = "Error bars = ±1 SE"
  )

ggsave("poster_fig1_rho.png", fig1, width = 7, height = 5, dpi = 300)
cat("Saved poster_fig1_rho.png\n")


source("R/simulation_source.R")
source("R/plot_traj.R")
source("R/plotting_source.R")

res <- readRDS("res.rds")
#fig1_blog <- 
  plot_metric_accuracy_v4(
  res,
  metric  = "rmse",
  family  = "correlated",
  xvar    = "rho",
  smooth  = TRUE,
  smooth_method = "loess",
  smooth_span   = 0.7,
  error_bar = FALSE,
  show_labels = FALSE
) +
  labs(title = NULL, subtitle = NULL,
       x = expression(Cross-country~alignment~(rho)),
       y = "Root Mean Square Error (RMSE)") +
  # Right-side labels
  annotate("text", x = 0.08, y = 0.355, label = "Fixed (20)",
           hjust = 0, size = 3.5, color = "#000000", family = "serif") +
  annotate("text", x = 0.08, y = 0.415, label = "CAT (rescoring)",
           hjust = 0, size = 3.5, color = "#E69F00", family = "serif") +
  annotate("text", x = 0.08, y = 0.315, label = "CAT (benchmark)",
           hjust = 0, size = 3.5, color = "#0072B2", family = "serif") +
  annotate("text", x = 0.08, y = 0.205, label = "Full (100)",
           hjust = 0, size = 3.5, color = "#CC79A7", family = "serif") +
  # Scoring error arrow (CAT-O to CAT-R at rho = 0.25)
  annotate("segment", x = 0.23, xend = 0.23,
           y = 0.63, yend = 0.38,
           arrow = arrow(ends = "both", length = unit(0.08, "inches")),
           color = "#D55E00", linewidth = 0.6) +
  annotate("text", x = 0.25, y = 0.51,
           label = "scoring error\n(fixable)", size = 5.2,
           color = "#D55E00", hjust = 0, family = "serif", fontface = "italic") +
  # Selection error arrow (CAT-R to CAT-B at rho = 0.25)
  annotate("segment", x = 0.23, xend = 0.23,
           y = 0.38, yend = 0.33,
           arrow = arrow(ends = "both", length = unit(0.08, "inches")),
           color = "#0072B2", linewidth = 0.6) +
  annotate("text", x = 0.25, y = 0.355,
           label = "selection error (not fixable)", size = 5.2,
           color = "#0072B2", hjust = 0, family = "serif", fontface = "italic") +
  # CAT-O label (high on the curve)
  annotate("text", x = 0.08, y = 0.72, label = "CAT (operational)",
           size = 3.5, color = "#D55E00", family = "serif") +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    expand = expansion(mult = c(0.02, 0.25))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    text             = element_text(family = "serif"),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#eeeeee", linewidth = 0.3),
    axis.title       = element_text(size = 13),
    axis.text        = element_text(size = 11, color = "#444444"),
    legend.position  = "none",
    strip.text       = element_blank()
  )

ggsave("blog_fig1_rho.png", fig1_blog,
       width = 8, height = 5, dpi = 300, bg = "white")
# ══════════════════════════════════════════════════════════════
# FIGURE 2: Stepwise RMSE trajectories at ρ = 0.75
# ══════════════════════════════════════════════════════════════
# Shows: CAT-O and CAT-R converge differently as test lengthens;
#        Fixed-form dots for comparison; Full-100 ceiling.
# n_reps=50 for smooth SE bands.

fig2 <- plot_traj_simple(
  res,
  bias_type    = "none",
  b_bias       = 0,
  err_mode     = "correlated",
  rho          = 0.75,
  n_persons    = 300,
  fixed_length = 50,
  max_steps    = 50,
  methods      = c("cat_ops", "cat_rescore", "full", "fixed_traj"),
  n_reps       = 50,
  seed         = 20250831,
  keep_only    = c("cat_ops", "cat_rescore", "fixed_traj", "full"),
  show_ribbon  = TRUE
)

ggsave("poster_fig2_stepwise.png", fig2, width = 7, height = 5, dpi = 300)
cat("Saved poster_fig2_stepwise.png\n")

cat("\nDone! Replace poster figures with these two PNGs.\n")

# ══════════════════════════════════════════════════════════════
# OPTIONAL: Additional figures for paper
# ══════════════════════════════════════════════════════════════

# # Additive family (σ): faceted by bias type
# fig_add <- plot_metric_accuracy_v4(
#   res, metric = "rmse", family = "additive", xvar = "sigma",
#   smooth = FALSE, error_bar = TRUE
# )
# ggsave("paper_fig_additive_rmse.png", fig_add, width = 12, height = 5, dpi = 300)

# # Pearson r vs ρ
# fig_r <- plot_metric_accuracy_v4(
#   res, metric = "pearson_r", family = "correlated", xvar = "rho",
#   smooth = TRUE
# )
# ggsave("paper_fig_correlated_r.png", fig_r, width = 7, height = 5, dpi = 300)

# # Stepwise at ρ = 0.50 (worse misalignment)
# fig2b <- plot_traj_simple(
#   res, bias_type = "none", err_mode = "correlated", rho = 0.50,
#   n_persons = 300, fixed_length = 50, max_steps = 50,
#   methods = c("cat_ops","cat_rescore","full","fixed_traj"),
#   n_reps = 50, seed = 20250831,
#   keep_only = c("cat_ops","cat_rescore","fixed_traj","full"),
#   show_ribbon = TRUE
# )
# ggsave("paper_fig_stepwise_rho50.png", fig2b, width = 7, height = 5, dpi = 300)
