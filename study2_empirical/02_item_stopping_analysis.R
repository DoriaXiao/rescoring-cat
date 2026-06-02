## ============================================================
## Item-level pattern + stopping-behavior analyses
##
## Addresses two questions Ben raised in the May 13 meeting:
##   Q1. Item-level pattern: did easy items (low b_old) have their
##       difficulties revised more than hard items? If so, this
##       explains why low-ability children shifted more.
##   Q2. Stopping behavior: do low-ability children stop sooner
##       in adaptive sessions? Links to the simulation's stepwise
##       prediction (rescoring gain largest at 5 to 15 items).
##
## USAGE: source("study2_empirical/01_rescoring_analysis.R") FIRST, then
## source this script in the SAME R session. It reuses these objects
## created by 01:
##   params      (from levante-runs/item parameter comparison.csv)
##   comparison  (per-child per-task with theta_old, theta_new_1pl,
##                age, site, run_id, n_items, diff_1pl)
##   run_data    (run-level metadata, includes adaptive flag)
##   table_dir, fig_dir
##
## Produces Appendix D outputs (item-level shift, stopping by ability).
## ============================================================

# ---- Task labels and colors (re-defined here so block is portable) ----
task_labels <- c(
  "egma-math"                = "Mathematics",
  "matrix-reasoning"         = "Matrix Reasoning",
  "mental-rotation"          = "Mental Rotation",
  "same-different-selection" = "Same-Different"
)
task_colors <- c(
  "Mathematics"      = "#2D5FA0",
  "Matrix Reasoning" = "#E8604C",
  "Mental Rotation"  = "#D4A843",
  "Same-Different"   = "#6A8E3C"
)


# ════════════════════════════════════════════════════════════
# Q1: ITEM-LEVEL PATTERN OF DIFFICULTY CHANGES
# ════════════════════════════════════════════════════════════
# Regress delta_b = b_new - b_old on b_old, per task.
# Negative slope: easy items revised harder more than hard items
# (Ben's hypothesis).
# Positive slope: opposite.

items_q1 <- params %>%
  filter(param_task_id %in% names(task_labels)) %>%
  mutate(delta_b   = b_new - b_old,
         task_name = recode(param_task_id, !!!task_labels))

cat("\n========== Q1: ITEM-LEVEL SHIFT ==========\n\n")
cat("Items per task with both b_old and b_new:\n")
count(items_q1, task_name) %>% print()

item_summary <- items_q1 %>%
  group_by(task_name) %>%
  summarise(
    n_items     = n(),
    r_old_new   = cor(b_old, b_new),
    mean_delta  = mean(delta_b),
    sd_delta    = sd(delta_b),
    slope_delta_on_bold = coef(lm(delta_b ~ b_old))[2],
    p_slope     = summary(lm(delta_b ~ b_old))$coefficients[2, "Pr(>|t|)"],
    .groups = "drop"
  )

cat("\n-- Item-level shift summary by task --\n")
print(item_summary)

write_csv(item_summary,
          file.path(table_dir, "results_item_level_shift.csv"))

# Plot: delta_b vs b_old, faceted by task
p_item <- ggplot(items_q1,
                 aes(x = b_old, y = delta_b)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
  geom_point(alpha = 0.4, size = 1.0, color = "grey30") +
  geom_smooth(method = "lm", se = TRUE,
              color = "black", fill = "grey70",
              linewidth = 0.6, alpha = 0.3) +
  facet_wrap(~ task_name, scales = "free") +
  labs(x = expression("Old difficulty ("*b[old]*")"),
       y = expression(Delta * b == b[new] - b[old])) +
  theme_bw(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        strip.background = element_rect(fill = "grey90", color = "grey50"),
        strip.text = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3))

ggsave(file.path(fig_dir, "fig_item_level_shift.png"), p_item,
       width = 7, height = 5.5, dpi = 300, bg = "white")
# ════════════════════════════════════════════════════════════
# Q2: STOPPING BEHAVIOR BY ABILITY TERCILE
# ════════════════════════════════════════════════════════════
# Use n_items directly from comparison (which inherits it from
# score_both, computed as trial count per scored run). Join only
# the adaptive flag from run_data.

comparison_with_adaptive <- comparison %>%
  left_join(run_data %>% select(run_id, adaptive),
            by = "run_id")

cat("\n========== Q2: STOPPING BEHAVIOR ==========\n\n")
cat("Adaptive vs non-adaptive by task:\n")
count(comparison_with_adaptive, task_id, adaptive) %>% print()

cat("\nn_items distribution (adaptive only):\n")
comparison_with_adaptive %>%
  filter(adaptive == TRUE) %>%
  group_by(task_id) %>%
  summarise(n          = n(),
            mean_items = mean(n_items, na.rm = TRUE),
            min_items  = min(n_items, na.rm = TRUE),
            max_items  = max(n_items, na.rm = TRUE),
            sd_items   = sd(n_items, na.rm = TRUE)) %>%
  print()

stop_data <- comparison_with_adaptive %>%
  filter(adaptive == TRUE, !is.na(n_items), !is.na(theta_new_1pl)) %>%
  group_by(task_id) %>%
  mutate(ability_group = factor(
    case_when(ntile(theta_new_1pl, 3) == 1 ~ "Low",
              ntile(theta_new_1pl, 3) == 2 ~ "Mid",
              ntile(theta_new_1pl, 3) == 3 ~ "High"),
    levels = c("Low", "Mid", "High"))) %>%
  ungroup() %>%
  mutate(task_name = recode(task_id, !!!task_labels))

stop_summary <- stop_data %>%
  group_by(task_name, ability_group) %>%
  summarise(n_children       = n(),
            mean_items       = mean(n_items),
            sd_items         = sd(n_items),
            mean_theta_shift = mean(diff_1pl),
            .groups = "drop")

cat("\n-- Items administered by ability tercile (adaptive only) --\n")
print(stop_summary, n = 30)

write_csv(stop_summary,
          file.path(table_dir, "results_stopping_by_ability.csv"))

# Per-task ANOVA
cat("\n-- ANOVA: n_items ~ ability_group --\n")
for (tk in unique(stop_data$task_id)) {
  d <- stop_data %>% filter(task_id == tk)
  if (n_distinct(d$ability_group) < 2 || nrow(d) < 6) {
    cat(sprintf("  %-26s skipped (n = %d, %d terciles)\n",
                recode(tk, !!!task_labels),
                nrow(d), n_distinct(d$ability_group)))
    next
  }
  fit <- aov(n_items ~ ability_group, data = d)
  cat(sprintf("  %-26s F = %5.2f, p = %.4f, n = %d\n",
              recode(tk, !!!task_labels),
              summary(fit)[[1]]["ability_group", "F value"],
              summary(fit)[[1]]["ability_group", "Pr(>F)"],
              nrow(d)))
}

# Plot
p_stop <- ggplot(stop_summary,
                 aes(x = ability_group, y = mean_items)) +
  geom_col(width = 0.6, fill = "grey50", color = "grey20",
           linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_items - sd_items / sqrt(n_children),
                    ymax = mean_items + sd_items / sqrt(n_children)),
                width = 0.18, color = "grey20", linewidth = 0.4) +
  facet_wrap(~ task_name, scales = "free_y") +
  labs(x = "Ability tercile (by rescored theta)",
       y = "Mean items administered") +
  theme_bw(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        strip.background = element_rect(fill = "grey90", color = "grey50"),
        strip.text = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3))

ggsave(file.path(fig_dir, "fig_stopping_by_ability.png"), p_stop,
       width = 7, height = 5.5, dpi = 300, bg = "white")

# ════════════════════════════════════════════════════════════
# CONNECTING Q1 AND Q2: does the item-level slope predict
# the Low-minus-High theta gap?
# ════════════════════════════════════════════════════════════

per_tercile_shift <- stop_data %>%
  group_by(task_name, ability_group) %>%
  summarise(mean_shift = mean(diff_1pl), .groups = "drop") %>%
  pivot_wider(names_from = ability_group, values_from = mean_shift) %>%
  mutate(low_minus_high = Low - High)

combined_q1_q2 <- item_summary %>%
  select(task_name, slope_delta_on_bold) %>%
  left_join(per_tercile_shift, by = "task_name")

cat("\n-- Q1 slope vs Q2 Low-minus-High shift --\n")
print(combined_q1_q2)

write_csv(combined_q1_q2,
          file.path(table_dir, "results_item_slope_vs_theta_gap.csv"))


cat("\n\nDone (Q1 + Q2).\n")
cat("Tables:\n")
cat("  ", file.path(table_dir, "results_item_level_shift.csv"), "\n")
cat("  ", file.path(table_dir, "results_stopping_by_ability.csv"), "\n")
cat("  ", file.path(table_dir, "results_item_slope_vs_theta_gap.csv"), "\n")
cat("Figures:\n")
cat("  ", file.path(fig_dir, "fig_item_level_shift.png"), "\n")
cat("  ", file.path(fig_dir, "fig_stopping_by_ability.png"), "\n")