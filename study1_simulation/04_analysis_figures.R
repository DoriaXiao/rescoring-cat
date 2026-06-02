## ============================================================
## Simulation Analysis v2: Figures & Tables from per-rep res.rds
##
## Reads:  output/data/res.rds        (per-rep, from v2 sim)
##         output/data/res_stepwise.rds (unchanged, from original)
##
## Outputs: same 7 figures + CSV tables as before, now with
##          SE bands on Figures 1, 2, 3.
##
## The per-rep res.rds has columns:
##   mode, rho, sigma, delta, condition, rep, rmse, bias, cor_r, mse
##
## This script aggregates to cell-level (mean +/- SE across reps)
## for plotting, and pools across all reps for paper-number checks.
## ============================================================

library(tidyverse)

res_raw  <- readRDS("output/data/res.rds")
stepwise <- readRDS("output/data/res_stepwise.rds")

table_dir <- "output/tables"
fig_dir   <- "figures/paper"
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,   showWarnings = FALSE, recursive = TRUE)

cat("res.rds:", nrow(res_raw), "rows (",
    n_distinct(res_raw$rep), "reps per cell )\n")
count(res_raw, mode) %>% print()
cat("stepwise:", nrow(stepwise), "rows\n")


# ══════════════════════════════════════════════════════════════
# 0. AGGREGATE: cell-level means + SE for plotting
# ══════════════════════════════════════════════════════════════

# Cell = unique combination of (mode, rho, sigma, delta, condition)
# Per-rep data has one row per (cell, rep)

res <- res_raw %>%
  group_by(mode, rho, sigma, delta, condition) %>%
  summarise(
    mean_rmse = mean(rmse),
    se_rmse   = sd(rmse) / sqrt(n()),
    mean_bias = mean(bias),
    mean_cor  = mean(cor_r),
    # Cell-level aggregates for paper numbers (pool across reps)
    cell_mse      = mean(mse),
    cell_rmse     = sqrt(mean(mse)),
    cell_bias     = mean(bias),
    cell_bias_sq  = mean(bias)^2,
    cell_pct_bias = if_else(mean(mse) > 1e-10,
                            mean(bias)^2 / mean(mse) * 100, 0),
    n_reps = n(),
    .groups = "drop"
  )


# ── Style constants ───────────────────────────────────────────

cond_order  <- c("CAT-O", "CAT-R", "CAT-B", "Fixed-20", "Full-100")
cond_colors <- c("CAT-O" = "#E8604C", "CAT-R" = "#D4A843",
                 "CAT-B" = "#2D5FA0", "Fixed-20" = "grey30",
                 "Full-100" = "#C48FB5")
cond_fills  <- cond_colors
cond_lines  <- c("CAT-O" = "solid", "CAT-R" = "solid",
                 "CAT-B" = "solid", "Fixed-20" = "dashed",
                 "Full-100" = "dotted")

base_theme <- theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.background = element_rect(fill = "white", color = NA))


# ══════════════════════════════════════════════════════════════
# 1. CROSS-POPULATION MISALIGNMENT (Figure 1)
# ══════════════════════════════════════════════════════════════

rho_data <- res %>% filter(mode == "correlated")

p_rmse_rho <- ggplot(rho_data, aes(x = rho, y = mean_rmse,
                                    color = condition, fill = condition,
                                    linetype = condition)) +
  geom_ribbon(aes(ymin = mean_rmse - se_rmse,
                  ymax = mean_rmse + se_rmse),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = cond_colors, breaks = cond_order) +
  scale_fill_manual(values = cond_fills, breaks = cond_order, guide = "none") +
  scale_linetype_manual(values = cond_lines, breaks = cond_order) +
  labs(x = expression("Cross-population alignment ("*rho*")"),
       y = "RMSE", color = NULL, linetype = NULL) +
  base_theme
ggsave(file.path(fig_dir, "fig_sim_rmse_rho.png"), p_rmse_rho,
       width = 7, height = 4.5, dpi = 300, bg = "white")

# Appendix: bias
p_bias_rho <- ggplot(rho_data, aes(x = rho, y = mean_bias,
                                    color = condition, linetype = condition)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
  scale_color_manual(values = cond_colors, breaks = cond_order) +
  scale_linetype_manual(values = cond_lines, breaks = cond_order) +
  labs(x = expression("Cross-population alignment ("*rho*")"),
       y = "Mean Bias", color = NULL, linetype = NULL) +
  base_theme
ggsave(file.path(fig_dir, "fig_sim_bias_rho.png"), p_bias_rho,
       width = 7, height = 4.5, dpi = 300, bg = "white")

# Appendix: correlation
p_cor_rho <- ggplot(rho_data, aes(x = rho, y = mean_cor,
                                   color = condition, linetype = condition)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1.2) +
  scale_color_manual(values = cond_colors, breaks = cond_order) +
  scale_linetype_manual(values = cond_lines, breaks = cond_order) +
  labs(x = expression("Cross-population alignment ("*rho*")"),
       y = expression("Correlation ("*r*")"), color = NULL, linetype = NULL) +
  base_theme
ggsave(file.path(fig_dir, "fig_sim_cor_rho.png"), p_cor_rho,
       width = 7, height = 4.5, dpi = 300, bg = "white")

write_csv(rho_data, file.path(table_dir, "results_sim_correlated.csv"))


# ══════════════════════════════════════════════════════════════
# 2. RANDOM NOISE (Figure 2)
# ══════════════════════════════════════════════════════════════

noise_data <- res %>% filter(mode == "random_noise")

p_rmse_noise <- ggplot(noise_data, aes(x = sigma, y = mean_rmse,
                                        color = condition, fill = condition,
                                        linetype = condition)) +
  geom_ribbon(aes(ymin = mean_rmse - se_rmse,
                  ymax = mean_rmse + se_rmse),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_color_manual(values = cond_colors, breaks = cond_order) +
  scale_fill_manual(values = cond_fills, breaks = cond_order, guide = "none") +
  scale_linetype_manual(values = cond_lines, breaks = cond_order) +
  labs(x = expression("Calibration noise ("*sigma*")"),
       y = "RMSE", color = NULL, linetype = NULL) +
  base_theme
ggsave(file.path(fig_dir, "fig_sim_rmse_noise.png"), p_rmse_noise,
       width = 7, height = 4.5, dpi = 300, bg = "white")

write_csv(noise_data, file.path(table_dir, "results_sim_noise.csv"))


# ══════════════════════════════════════════════════════════════
# 3. STRUCTURED SHIFTS (Figure 3, 3-panel)
# ══════════════════════════════════════════════════════════════

shift_data <- res %>% filter(mode == "structured_shift")
delta_labels <- c("-0.3" = "\u0394b = \u22120.3",
                  "0.3"  = "\u0394b = +0.3",
                  "1"    = "\u0394b = +1.0")

p_rmse_shift <- ggplot(shift_data, aes(x = sigma, y = mean_rmse,
                                        color = condition, fill = condition,
                                        linetype = condition)) +
  geom_ribbon(aes(ymin = mean_rmse - se_rmse,
                  ymax = mean_rmse + se_rmse),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  facet_wrap(~delta, labeller = labeller(delta = delta_labels)) +
  scale_color_manual(values = cond_colors, breaks = cond_order) +
  scale_fill_manual(values = cond_fills, breaks = cond_order, guide = "none") +
  scale_linetype_manual(values = cond_lines, breaks = cond_order) +
  labs(x = expression("Calibration noise ("*sigma*")"),
       y = "RMSE", color = NULL, linetype = NULL) +
  base_theme +
  theme(strip.text = element_text(face = "bold"))
ggsave(file.path(fig_dir, "fig_sim_rmse_shift.png"), p_rmse_shift,
       width = 10, height = 4.5, dpi = 300, bg = "white")

write_csv(shift_data, file.path(table_dir, "results_sim_shift.csv"))


# ══════════════════════════════════════════════════════════════
# 4. STEPWISE TRAJECTORY (Figure 4, unchanged data)
# ══════════════════════════════════════════════════════════════

full100_rmse <- res %>%
  filter(mode == "correlated", condition == "Full-100") %>%
  #filter(rho > 0.7 & rho < 0.8) %>%
  pull(cell_rmse) %>% mean()

step_colors <- c("CAT-O" = "#E8604C", "CAT-R" = "#D4A843",
                 "CAT-B" = "#2D5FA0", "Fixed-20" = "grey30")

p_stepwise <- ggplot(stepwise, aes(x = step, y = rmse,
                                    color = condition)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = full100_rmse, linetype = "dotted",
             color = "#C48FB5", linewidth = 0.8) +
  annotate("text", x = 2, y = full100_rmse - 0.02, label = "Full-100",
           color = "#C48FB5", size = 3.5, hjust = 0) +
  scale_color_manual(values = step_colors,
                     breaks = c("CAT-O", "CAT-R", "CAT-B", "Fixed-20")) +
  labs(x = "Number of items administered",
       y = "RMSE", color = NULL,
       caption = expression("Stepwise RMSE at "*rho*" = 0.75, n = 300, 50 replications")) +
  base_theme
ggsave(file.path(fig_dir, "fig_sim_stepwise.png"), p_stepwise,
       width = 7, height = 4.5, dpi = 300, bg = "white")

write_csv(stepwise, file.path(table_dir, "results_sim_stepwise.csv"))


# ══════════════════════════════════════════════════════════════
# 5. PAPER NUMBER VERIFICATION
#    Print each number that appears in the manuscript text,
#    with its section reference, so you can update if changed.
# ══════════════════════════════════════════════════════════════

cat("\n\n")
cat("══════════════════════════════════════════════════════════\n")
cat("  PAPER NUMBER VERIFICATION\n")
cat("  Compare these to your manuscript and update as needed.\n")
cat("══════════════════════════════════════════════════════════\n\n")

# --- 5a. Scoring error as % of total error ---
# Paper: "76% (cross-pop), 85% (random), 97% (structured)"
# Locations: Abstract, Section 2.4.5, Section 4.1

scoring_pct <- res %>%
  select(mode, rho, sigma, delta, condition, cell_rmse) %>%
  distinct() %>%
  pivot_wider(names_from = condition, values_from = cell_rmse) %>%
  mutate(
    total_error   = `CAT-O` - `CAT-B`,
    scoring_error = `CAT-O` - `CAT-R`,
    pct_scoring   = if_else(total_error > 1e-6,
                            scoring_error / total_error * 100, NA_real_)
  ) %>%
  group_by(mode) %>%
  summarise(
    mean_pct = mean(pct_scoring, na.rm = TRUE),
    min_pct  = min(pct_scoring, na.rm = TRUE),
    .groups  = "drop"
  )

cat("--- Scoring error % of total (Abstract, Sec 2.4.5, Sec 4.1) ---\n")
cat("  Paper says: 76% (cross-pop), 85% (random), 97% (structured)\n")
cat("  New values:\n")
scoring_pct %>% print()


# --- 5b. Bias decomposition (Table 2) ---
# Paper: "bias < 0.5% of MSE after rescoring across all error families"
# Location: Abstract, Section 2.4.3, Table 2, Section 4.1

bias_decomp <- res %>%
  select(mode, rho, sigma, delta, condition, cell_pct_bias) %>%
  distinct() %>%
  filter(condition %in% c("CAT-O", "CAT-R")) %>%
  group_by(mode, condition) %>%
  summarise(
    mean_pct_bias = mean(cell_pct_bias),
    max_pct_bias  = max(cell_pct_bias),
    .groups = "drop"
  )

cat("\n--- Bias decomposition for Table 2 ---\n")
cat("  Paper says: CAT-R max < 0.5% across all families\n")
cat("  New values:\n")
bias_decomp %>%
  pivot_wider(names_from = condition,
              values_from = c(mean_pct_bias, max_pct_bias)) %>%
  print()

max_catr_bias <- max(bias_decomp$max_pct_bias[bias_decomp$condition == "CAT-R"])
cat(sprintf("\n  Max CAT-R pct_bias: %.4f%% (paper claims < 0.5%%: %s)\n",
            max_catr_bias,
            if_else(max_catr_bias < 0.5, "CONFIRMED", "CHECK")))


# --- 5c. Structured shift CAT-O max bias ---
# Paper: "up to 90% of MSE at delta = 1.0"
# Location: Section 2.4.3

max_shift_bias <- res %>%
  select(mode, delta, condition, cell_pct_bias) %>%
  distinct() %>%
  filter(mode == "structured_shift", condition == "CAT-O",
         abs(delta - 1.0) < 0.01) %>%
  summarise(max_pct = max(cell_pct_bias)) %>%
  pull(max_pct)

cat(sprintf("\n--- Structured shift CAT-O max bias (Sec 2.4.3) ---\n"))
cat(sprintf("  Paper says: up to 90%% of MSE at delta = 1.0\n"))
cat(sprintf("  New value: %.1f%%\n", max_shift_bias))


# --- 5d. Cross-population & random noise bias "typically < 5%" ---
# Location: Section 2.4.3

cat("\n--- CAT-O bias 'typically < 5%' for correlated & noise (Sec 2.4.3) ---\n")
res %>%
  select(mode, condition, cell_pct_bias) %>%
  distinct() %>%
  filter(condition == "CAT-O", mode %in% c("correlated", "random_noise")) %>%
  group_by(mode) %>%
  summarise(mean_pct = mean(cell_pct_bias),
            max_pct = max(cell_pct_bias),
            .groups = "drop") %>%
  print()


# --- 5e. Full Table 2 for the paper ---
cat("\n--- FULL TABLE 2 (copy into paper) ---\n")
bias_decomp %>%
  mutate(
    cell = paste0(sprintf("%.1f%%", mean_pct_bias), " / ",
                  sprintf("%.1f%%", max_pct_bias))
  ) %>%
  select(mode, condition, cell) %>%
  pivot_wider(names_from = condition, values_from = cell) %>%
  print()


# ══════════════════════════════════════════════════════════════
# 6. SAVE CELL-LEVEL TABLE
# ══════════════════════════════════════════════════════════════

write_csv(
  res %>% select(mode, rho, sigma, delta, condition,
                 cell_rmse, cell_bias, cell_mse, cell_pct_bias) %>%
    distinct(),
  file.path(table_dir, "results_sim_bias_decomposition.csv")
)

write_csv(
  scoring_pct,
  file.path(table_dir, "results_sim_scoring_error_pct.csv")
)


# ══════════════════════════════════════════════════════════════
# DONE
# ══════════════════════════════════════════════════════════════

cat("\n\nDone!\n")
cat("Figures:\n")
cat(paste(" ", list.files(fig_dir, "\\.png$")), sep = "\n")
cat("\nTables:\n")
cat(paste(" ", list.files(table_dir, "results_sim")), sep = "\n")

# ══════════════════════════════════════════════════════════════
# 1. CROSS-POPULATION MISALIGNMENT (Figure 1)
# ══════════════════════════════════════════════════════════════

rho_data <- res %>% filter(mode == "correlated")

# Anchor the scoring/selection error annotations at a rho level where
# the CAT-O vs CAT-R gap is clearly visible. We pick the rho closest to
# 0.22 (close to the blog figure's anchor point).
anchor_rho <- rho_data$rho[which.min(abs(rho_data$rho - 0.22))]

anchor_vals <- rho_data %>%
  filter(abs(rho - anchor_rho) < 1e-6,
         condition %in% c("CAT-O", "CAT-R", "CAT-B")) %>%
  select(condition, mean_rmse) %>%
  tibble::deframe()

# Midpoints for label placement
scoring_mid   <- (anchor_vals[["CAT-O"]] + anchor_vals[["CAT-R"]]) / 2
selection_mid <- (anchor_vals[["CAT-R"]] + anchor_vals[["CAT-B"]]) / 2

p_rmse_rho <- 
  ggplot(rho_data, aes(x = rho, y = mean_rmse,
                                   color = condition, fill = condition,
                                   linetype = condition)) +
  geom_ribbon(aes(ymin = mean_rmse - se_rmse,
                  ymax = mean_rmse + se_rmse),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  # Vertical braces marking the two error components at anchor_rho
  annotate("segment",
           x = anchor_rho, xend = anchor_rho,
           y = anchor_vals[["CAT-R"]], yend = anchor_vals[["CAT-O"]],
           arrow = arrow(ends = "both", length = unit(0.15, "cm")),
           color = "grey30", linewidth = 0.4) +
  annotate("segment",
           x = anchor_rho, xend = anchor_rho,
           y = anchor_vals[["CAT-B"]], yend = anchor_vals[["CAT-R"]],
           arrow = arrow(ends = "both", length = unit(0.15, "cm")),
           color = "grey30", linewidth = 0.4) +
  # Italic labels offset to the right of the braces
  annotate("text",
           x = anchor_rho - 0.1, y = scoring_mid,
           label = "scoring error\n(fixable)",
           hjust = 0, size = 4, fontface = "italic",
           color = "grey20", lineheight = 0.9) +
  annotate("text",
           x = anchor_rho + 0.03, y = selection_mid-0.01,
           label = "selection error (not fixable)",
           hjust = 0, size = 4, fontface = "italic",
           color = "grey20", lineheight = 0.9) +
  scale_color_manual(values = cond_colors, breaks = cond_order) +
  scale_fill_manual(values = cond_fills, breaks = cond_order, guide = "none") +
  scale_linetype_manual(values = cond_lines, breaks = cond_order) +
  labs(x = expression("Cross-population alignment ("*rho*")"),
       y = "RMSE", color = NULL, linetype = NULL) +
  base_theme
ggsave(file.path(fig_dir, "fig_sim_rmse_rho.png"), p_rmse_rho,
       width = 7, height = 4.5, dpi = 300, bg = "white")

## ============================================================
## PATCH: Figure 1 — replace arrows with shaded regions
## ============================================================
## Per Mike's feedback ("should these be shaded regions, not
## distances?"). Shaded regions correctly show error components
## across the full rho range; arrows only show the gap at one
## anchor point.
##
## WHERE TO PUT THIS:
##   Replace the SECOND "1. CROSS-POPULATION MISALIGNMENT" block
##   in run_simulation_analysis.R (the one with anchor_rho and
##   the annotate("segment", ...) arrows). This block overwrites
##   fig_sim_rmse_rho.png, so we want it to be the final version.
##
## ASSUMES (already defined earlier in your script):
##   res           — aggregated tibble with mean_rmse, se_rmse
##   cond_colors   — named vector with CAT-O, CAT-R, CAT-B, etc.
##   cond_fills, cond_lines, cond_order
##   base_theme
##   fig_dir, table_dir
## ============================================================

# ══════════════════════════════════════════════════════════════
# 1. CROSS-POPULATION MISALIGNMENT (Figure 1, shaded regions)
# ══════════════════════════════════════════════════════════════
rho_data <- res %>% filter(mode == "correlated")
# ---- Pivot wide to build ribbon bounds ----
rho_wide <- rho_data %>%
  select(rho, condition, mean_rmse) %>%
  pivot_wider(names_from = condition, values_from = mean_rmse)
stopifnot(all(c("CAT-O", "CAT-R", "CAT-B") %in% names(rho_wide)))
# Scoring error region: between CAT-O (upper) and CAT-R (lower)
ribbon_scoring <- rho_wide %>%
  transmute(rho, ymin = `CAT-R`, ymax = `CAT-O`)
# Selection error region: between CAT-R (upper) and CAT-B (lower)
ribbon_selection <- rho_wide %>%
  transmute(rho, ymin = `CAT-B`, ymax = `CAT-R`)
# ---- Anchor points for the region labels ----
# Scoring label at rho ≈ 0.22 (left where the gap is widest)
# Selection label at rho ≈ 0.55 (middle, comfortable spot)
label_anchor_scoring   <- rho_wide$rho[which.min(abs(rho_wide$rho - 0.22))]
label_anchor_selection <- rho_wide$rho[which.min(abs(rho_wide$rho - 0.55))]
scoring_y   <- with(rho_wide %>% filter(rho == label_anchor_scoring),
                    mean(c(`CAT-O`, `CAT-R`)))
selection_y <- with(rho_wide %>% filter(rho == label_anchor_selection),
                    mean(c(`CAT-R`, `CAT-B`)))
# ---- Build the figure ----
p_rmse_rho <- ggplot(rho_data, aes(x = rho, y = mean_rmse,
                                   color = condition, fill = condition,
                                   linetype = condition)) +
  # Error-component shaded regions (drawn first, underneath lines)
  # Convention: each region inherits the color of its upper-bounding line.
  # Scoring region (CAT-O to CAT-R): filled red, matching CAT-O.
  # Selection region (CAT-R to CAT-B): filled gold, matching CAT-R.
  geom_ribbon(data = ribbon_scoring,
              inherit.aes = FALSE,
              aes(x = rho, ymin = ymin, ymax = ymax),
              fill = cond_colors[["CAT-O"]], alpha = 0.18) +
  geom_ribbon(data = ribbon_selection,
              inherit.aes = FALSE,
              aes(x = rho, ymin = ymin, ymax = ymax),
              fill = cond_colors[["CAT-R"]], alpha = 0.30) +
  
  # Lines and points (on top of ribbons)
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  
  # Italic labels inside each region
  annotate("text",
           x = label_anchor_scoring, y = scoring_y,
           label = "scoring error\n(fixable)",
           hjust = 0.5, size = 4, fontface = "italic",
           color = cond_colors[["CAT-O"]], lineheight = 0.9) +
  annotate("text",
           x = label_anchor_selection, y = selection_y,
           label = "selection error (not fixable)",
           hjust = 0.5, size = 3.6, fontface = "italic",
           color = cond_colors[["CAT-R"]], lineheight = 0.9) +
  
  scale_color_manual(name = NULL,
                     values = cond_colors,
                     breaks = cond_order,
                     labels = cond_order) +
  scale_fill_manual(name = NULL,
                    values = cond_fills,
                    breaks = cond_order,
                    labels = cond_order,
                    guide = "none") +
  scale_linetype_manual(name = NULL,
                        values = cond_lines,
                        breaks = cond_order,
                        labels = cond_order) +
  labs(x = expression("Cross-population alignment (" * rho * ")"),
       y = "RMSE") +
  base_theme
p_rmse_rho <- p_rmse_rho +
  guides(color    = guide_legend(title = NULL, override.aes = list(
    linetype = cond_lines[cond_order],
    shape    = c(16, 16, 16, NA, NA))),
    linetype = "none",
    fill     = "none")
ggsave(file.path(fig_dir, "fig_sim_rmse_rho.png"), p_rmse_rho,
       width = 7, height = 4.5, dpi = 300, bg = "white")
