## ============================================================
## Plot perturbed-parameter sensitivity sim (Appendix D)
##
## Input:
##   output/data/res_perturbed.rds (per-examinee, from run_perturbed_sim.R)
##
## Output:
##   figures/paper/fig_sim_perturbed.png
##   output/tables/results_sim_perturbed.csv
##
## STYLE: Matches the 3-panel structured-shift figure in
## run_simulation_analysis_v2.R (lines 161-178):
##   - width = 10, height = 4.5, dpi = 300, bg = "white"
##   - shared cond_colors / cond_lines palette
##   - geom_ribbon for SE bands (alpha = 0.15)
##   - geom_line linewidth = 1, geom_point size = 1.5
##   - facet_wrap with bold strip text
##   - NO caption in figure (caption goes in LaTeX)
##
## CAT-Rhat is added as a new condition with a darker gold color and
## a dashed line, distinguishable from CAT-R (gold solid) and Fixed-20
## (grey dashed).
##
## Usage: run from the repo root (open rescoring-cat.Rproj), then source this script.
## Runtime: < 5 seconds.
## ============================================================

library(tidyverse)

# ── 1. Read data ──────────────────────────────────────────────

res_perturbed <- readRDS("output/data/res_perturbed.rds")
cat("Loaded res_perturbed.rds:", nrow(res_perturbed), "rows\n")

# ── 2. Aggregate to cell level (per-rep RMSE -> mean ± SE) ────

per_rep <- res_perturbed %>%
  group_by(rho, N_rescore, condition, rep) %>%
  summarise(rmse_rep = sqrt(mean((theta_hat - theta_true)^2)),
            bias_rep = mean(theta_hat - theta_true),
            .groups = "drop")

cell <- per_rep %>%
  group_by(rho, N_rescore, condition) %>%
  summarise(mean_rmse = mean(rmse_rep),
            se_rmse   = sd(rmse_rep) / sqrt(n()),
            mean_bias = mean(bias_rep),
            .groups = "drop")

# ── 3. Style — extend the existing palette with CAT-Rhat ─────

# Match cond_colors / cond_lines from run_simulation_analysis_v2.R lines 62-68,
# extended with a 5th condition CAT-Rhat. Color is a darker gold so it sits
# next to CAT-R (#D4A843) without clashing.
cond_colors <- c("CAT-O"    = "#E8604C",
                 "CAT-R"    = "#D4A843",
                 "CAT-Rhat" = "#8B6914",  # darker gold
                 "CAT-B"    = "#2D5FA0")
cond_lines  <- c("CAT-O"    = "solid",
                 "CAT-R"    = "solid",
                 "CAT-Rhat" = "dashed",
                 "CAT-B"    = "solid")
cond_order  <- c("CAT-O", "CAT-R", "CAT-Rhat", "CAT-B")

base_theme <- theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.background = element_rect(fill = "white", color = NA))

# ── 4. Plot ───────────────────────────────────────────────────

# Set condition factor order (matches main figures' legend ordering)
cell$condition <- factor(cell$condition, levels = cond_order)

# Facet labels parsed as math expression (rho = 0.34 etc.)
cell <- cell %>%
  mutate(rho_label = factor(rho,
                             levels = c(0.34, 0.50, 0.75),
                             labels = c("rho == 0.34",
                                        "rho == 0.50",
                                        "rho == 0.75")))

p_perturbed <- ggplot(cell, aes(x = N_rescore, y = mean_rmse,
                                 color = condition, fill = condition,
                                 linetype = condition)) +
  geom_ribbon(aes(ymin = mean_rmse - se_rmse,
                  ymax = mean_rmse + se_rmse),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~ rho_label, labeller = label_parsed) +
  scale_x_log10(breaks = c(100, 300, 1000)) +
  scale_color_manual(values = cond_colors, breaks = cond_order) +
  scale_fill_manual(values = cond_colors, breaks = cond_order, guide = "none") +
  scale_linetype_manual(values = cond_lines, breaks = cond_order) +
  labs(x = expression("Calibration sample size ("*italic(N)[rescore]*")"),
       y = "RMSE", color = NULL, linetype = NULL) +
  base_theme +
  theme(strip.text = element_text(face = "bold"))

# ── 5. Save ───────────────────────────────────────────────────

fig_dir   <- "figures/paper"
table_dir <- "output/tables"
dir.create(fig_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(fig_dir, "fig_sim_perturbed.png"), p_perturbed,
       width = 10, height = 4.5, dpi = 300, bg = "white")
write_csv(cell, file.path(table_dir, "results_sim_perturbed.csv"))

cat("\nWrote:\n")
cat("  ", file.path(fig_dir, "fig_sim_perturbed.png"), "\n")
cat("  ", file.path(table_dir, "results_sim_perturbed.csv"), "\n")

# ── 6. Headline numbers (for Appendix D paragraph) ────────────

cat("\n=== Headline numbers ===\n")
diff_summary <- cell %>%
  select(rho, N_rescore, condition, mean_rmse) %>%
  pivot_wider(names_from = condition, values_from = mean_rmse) %>%
  mutate(rhat_minus_r = `CAT-Rhat` - `CAT-R`,
         o_minus_rhat = `CAT-O` - `CAT-Rhat`)

cat(sprintf("Max |CAT-Rhat - CAT-R| RMSE: %.4f logits\n",
            max(abs(diff_summary$rhat_minus_r))))
cat(sprintf("Min RMSE reduction CAT-O -> CAT-Rhat: %.3f logits\n",
            min(diff_summary$o_minus_rhat)))
cat(sprintf("Max RMSE reduction CAT-O -> CAT-Rhat: %.3f logits\n",
            max(diff_summary$o_minus_rhat)))

bias_summ <- res_perturbed %>%
  group_by(rho, N_rescore, condition) %>%
  summarise(bias = mean(theta_hat - theta_true),
            mse  = mean((theta_hat - theta_true)^2),
            bias_pct_mse = bias^2 / mse * 100,
            .groups = "drop")

cat(sprintf("Max CAT-Rhat bias %% of MSE: %.3f %%\n",
            max(bias_summ$bias_pct_mse[bias_summ$condition == "CAT-Rhat"])))
