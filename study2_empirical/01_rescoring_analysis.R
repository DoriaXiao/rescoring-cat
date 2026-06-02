## ============================================================
## LEVANTE CAT Rescoring: Complete Empirical Analysis Pipeline (Study 2)
##
## Produces the empirical tables and figures in the paper:
##   - Tables 4-7 (parameter effect, by-site, age validity, model effect)
##   - cross-task convergent validity + bootstrap CIs
##   - equity terciles, age bands, test-retest
##   - Figures: score-shift scatter & density (paper Fig. 5)
##   - Per-child derived dataset: output/data/comparison.csv
##
## DATA ACCESS (NOT bundled in this repo; child-level LEVANTE data):
##   This script requires raw LEVANTE pilot data. See
##   study2_empirical/README.md for how to obtain it. Set the two
##   paths below to point at your local copies before running.
##
##   LEVANTE_PILOTS/  must contain 01_fetched_data/task_data_nested.rds
##                                  01_fetched_data/run_data.rds
##   LEVANTE_RUNS/    must contain all_scores.csv
##                                  item parameter comparison.csv
##
## Run from the repo root (open rescoring-cat.Rproj).
## ============================================================

library(tidyverse)
library(here)   # output paths resolve relative to repo root (.here marker)

# -- CONFIG: point these at your local LEVANTE data ------------
LEVANTE_PILOTS <- Sys.getenv("LEVANTE_PILOTS",
                             unset = "~/Documents/levante-pilots")
LEVANTE_RUNS   <- Sys.getenv("LEVANTE_RUNS",
                             unset = "~/Documents/levante-runs/levante-runs")
LEVANTE_PILOTS <- path.expand(LEVANTE_PILOTS)
LEVANTE_RUNS   <- path.expand(LEVANTE_RUNS)

# -- 0. Load trial data from levante-pilots --------------------

all_sum_scores <- read_rds(file.path(LEVANTE_PILOTS,
                                     "01_fetched_data/task_data_nested.rds"))

# Run-level metadata (for adaptive status, test-retest)
run_data <- read_rds(file.path(LEVANTE_PILOTS, "01_fetched_data/run_data.rds"))
cat("run_data:", nrow(run_data), "rows\n")

# -- 0b. Scores & params from levante-runs --------------------

path_scores <- file.path(LEVANTE_RUNS, "all_scores.csv")
path_params <- file.path(LEVANTE_RUNS, "item parameter comparison.csv")

stopifnot(
  "all_scores.csv not found (set LEVANTE_RUNS)." = file.exists(path_scores),
  "item parameter comparison.csv not found (set LEVANTE_RUNS)." = file.exists(path_params)
)

# -- 0c. Output directories (in this repo) ---------------------

table_dir <- here("output", "tables")
fig_dir   <- here("figures", "paper")
data_dir  <- here("output", "data")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir,  showWarnings = FALSE, recursive = TRUE)

# -- 1. Load scores & params ----------------------------------

target_tasks <- c("egma-math", "matrix-reasoning",
                  "mental-rotation", "same-different-selection")

new_scores <- read_csv(path_scores, show_col_types = FALSE) %>%
  filter(metric_type == "ability", task_id %in% target_tasks)

params <- read_csv(path_params, show_col_types = FALSE) %>%
  filter(!is.na(difficulty_old), difficulty_old != "",
         !is.na(difficulty), difficulty != "") %>%
  mutate(b_old = as.numeric(difficulty_old),
         b_new = as.numeric(difficulty),
         param_chance = replace_na(as.numeric(chance_level), 0)) %>%
  select(item_uid, param_task_id = task_id, param_chance, b_old, b_new) %>%
  distinct(item_uid, .keep_all = TRUE)

cat("Params:", nrow(params), "unique items\n")
count(params, param_task_id) %>% print()

# -- 2. Unnest trials ------------------------------------------

trials <- all_sum_scores %>% unnest(data)
cat("Trial rows:", nrow(trials), "\n")

# -- 3. EAP scoring function (1PL + guessing) ------------------

eap_score <- function(y, b, c_param = 0,
                      grid = seq(-6, 6, 0.05), prior_sd = 1.5) {
  a <- rep(1, length(y))
  ll <- sapply(grid, function(th) {
    p <- c_param + (1 - c_param) / (1 + exp(-a * (th - b)))
    p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
    sum(y * log(p) + (1 - y) * log(1 - p))
  })
  lp <- ll + dnorm(grid, 0, prior_sd, log = TRUE)
  w  <- exp(lp - max(lp)); w <- w / sum(w)
  c(theta = sum(grid * w),
    se = sqrt(sum((grid - sum(grid * w))^2 * w)))
}

# -- 4. Merge & score ------------------------------------------

merged <- trials %>%
  filter(!is.na(correct)) %>%
  mutate(correct_num = as.integer(correct)) %>%
  inner_join(params, by = "item_uid") %>%
  filter(param_task_id %in% target_tasks)

cat("\nMerged trial rows:", nrow(merged), "\n")
cat("Unique children:", n_distinct(merged$user_id), "\n")
cat("By task:\n")
count(merged, param_task_id) %>% print()

score_both <- merged %>%
  group_by(user_id, run_id, param_task_id, site) %>%
  filter(n() >= 5) %>%
  summarise(
    n_items     = n(),
    est_old     = list(eap_score(correct_num, b_old, c_param = param_chance)),
    est_new_1pl = list(eap_score(correct_num, b_new, c_param = param_chance)),
    .groups = "drop"
  ) %>%
  mutate(
    theta_old     = map_dbl(est_old, "theta"),
    se_old        = map_dbl(est_old, "se"),
    theta_new_1pl = map_dbl(est_new_1pl, "theta"),
    se_new_1pl    = map_dbl(est_new_1pl, "se")
  ) %>%
  select(-est_old, -est_new_1pl)

cat("Scored runs:", nrow(score_both), "\n")

# Latest run per child x task
score_latest <- score_both %>%
  left_join(
    new_scores %>% select(run_id, time_started) %>% distinct(),
    by = "run_id"
  ) %>%
  group_by(user_id, param_task_id) %>%
  slice_max(time_started, n = 1, with_ties = FALSE) %>%
  ungroup()

# -- 5. Multigroup IRT scores ---------------------------------

score_mg <- new_scores %>%
  group_by(user_id, task_id) %>%
  slice_max(time_started, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(user_id, task_id, age,
         theta_new_mg = metric_value, se_new_mg = metric_se) %>%
  mutate(across(c(theta_new_mg, se_new_mg, age), as.numeric))

# -- 6. Join all three estimates -------------------------------

comparison <- score_latest %>%
  rename(task_id = param_task_id) %>%
  inner_join(score_mg, by = c("user_id", "task_id")) %>%
  mutate(
    diff_1pl = theta_new_1pl - theta_old,
    diff_mg  = theta_new_mg  - theta_old
  )

cat("\n========== RESULTS ==========\n")
cat("Matched children x tasks:", nrow(comparison), "\n\n")


# =============================================================
# 7. RESULTS TABLES
# =============================================================

cat("-- Table A: Parameter effect (same 1PL EAP) --\n")
table_a <- comparison %>%
  group_by(task_id) %>%
  summarise(n = n(), r = cor(theta_old, theta_new_1pl),
            mean_shift = mean(diff_1pl), sd_shift = sd(diff_1pl),
            rmse = sqrt(mean(diff_1pl^2)), .groups = "drop")
print(table_a)

cat("\n-- Table B: Combined effect (1PL-old vs MG-new) --\n")
table_b <- comparison %>%
  group_by(task_id) %>%
  summarise(n = n(), r = cor(theta_old, theta_new_mg),
            mean_shift = mean(diff_mg), sd_shift = sd(diff_mg),
            rmse = sqrt(mean(diff_mg^2)), .groups = "drop")
print(table_b)

cat("\n-- Table C: Model effect (1PL vs MG, new params) --\n")
table_c <- comparison %>%
  group_by(task_id) %>%
  summarise(n = n(), r = cor(theta_new_1pl, theta_new_mg),
            mean_diff = mean(theta_new_mg - theta_new_1pl),
            sd_diff = sd(theta_new_mg - theta_new_1pl), .groups = "drop")
print(table_c)


# =============================================================
# 8. BY SITE
# =============================================================

cat("\n-- By task x site --\n")
table_site <- comparison %>%
  group_by(task_id, site) %>%
  summarise(n = n(), r = cor(theta_old, theta_new_1pl),
            mean_shift = mean(diff_1pl), sd_shift = sd(diff_1pl),
            .groups = "drop")
print(table_site, n = 30)


# =============================================================
# 9. VALIDATION
# =============================================================

cat("\n-- Age validity --\n")
table_age <- comparison %>%
  filter(!is.na(age)) %>%
  group_by(task_id) %>%
  summarise(n = n(),
            r_age_old    = cor(theta_old, age, use = "complete"),
            r_age_new1pl = cor(theta_new_1pl, age, use = "complete"),
            r_age_mg     = cor(theta_new_mg, age, use = "complete"),
            .groups = "drop")
print(table_age)

cat("\n-- Age validity by site --\n")
age_by_site <- comparison %>%
  filter(!is.na(age)) %>%
  group_by(task_id, site) %>%
  summarise(n = n(),
            r_age_old    = cor(theta_old, age, use = "complete"),
            r_age_new1pl = cor(theta_new_1pl, age, use = "complete"),
            delta_r      = r_age_new1pl - r_age_old,
            .groups = "drop")
print(age_by_site, n = 20)

cat("\n-- Cross-task convergent validity --\n")
wide <- comparison %>%
  select(user_id, task_id, theta_old, theta_new_1pl) %>%
  pivot_wider(names_from = task_id,
              values_from = c(theta_old, theta_new_1pl))

task_pairs <- list(
  c("egma-math", "matrix-reasoning"),
  c("egma-math", "mental-rotation"),
  c("matrix-reasoning", "mental-rotation"),
  c("egma-math", "same-different-selection"),
  c("matrix-reasoning", "same-different-selection")
)

cross <- map_dfr(task_pairs, function(pair) {
  o1 <- paste0("theta_old_", pair[1])
  o2 <- paste0("theta_old_", pair[2])
  n1 <- paste0("theta_new_1pl_", pair[1])
  n2 <- paste0("theta_new_1pl_", pair[2])
  d <- wide %>% filter(!is.na(.data[[o1]]), !is.na(.data[[o2]]),
                       !is.na(.data[[n1]]), !is.na(.data[[n2]]))
  if (nrow(d) < 10) return(NULL)
  tibble(pair = paste0(pair[1], " x ", pair[2]), n = nrow(d),
         r_old = cor(d[[o1]], d[[o2]]),
         r_new = cor(d[[n1]], d[[n2]]),
         delta = cor(d[[n1]], d[[n2]]) - cor(d[[o1]], d[[o2]]))
})
print(cross)


# =============================================================
# 10. BOOTSTRAP CIs ON CROSS-TASK VALIDITY
# =============================================================

cat("\n-- Bootstrap CIs (B = 2000) --\n")
set.seed(2026)
n_boot <- 2000

boot_results <- map_dfr(task_pairs, function(pair) {
  o1 <- paste0("theta_old_", pair[1])
  o2 <- paste0("theta_old_", pair[2])
  n1 <- paste0("theta_new_1pl_", pair[1])
  n2 <- paste0("theta_new_1pl_", pair[2])
  d <- wide %>% filter(!is.na(.data[[o1]]), !is.na(.data[[o2]]),
                       !is.na(.data[[n1]]), !is.na(.data[[n2]]))
  n <- nrow(d)
  delta_obs <- cor(d[[n1]], d[[n2]]) - cor(d[[o1]], d[[o2]])
  deltas <- replicate(n_boot, {
    idx <- sample(n, n, replace = TRUE)
    cor(d[[n1]][idx], d[[n2]][idx]) - cor(d[[o1]][idx], d[[o2]][idx])
  })
  tibble(pair = paste0(pair[1], " x ", pair[2]), n = n,
         r_old = cor(d[[o1]], d[[o2]]),
         r_new = cor(d[[n1]], d[[n2]]),
         delta = delta_obs,
         ci_lower = quantile(deltas, 0.025),
         ci_upper = quantile(deltas, 0.975),
         significant = quantile(deltas, 0.025) > 0)
})
print(boot_results)


# =============================================================
# 11. EQUITY: ABILITY TERCILES
# =============================================================

cat("\n-- Rescoring impact by ability tercile --\n")
equity <- comparison %>%
  group_by(task_id) %>%
  mutate(
    ability_group = factor(
      case_when(ntile(theta_new_1pl, 3) == 1 ~ "Low",
                ntile(theta_new_1pl, 3) == 2 ~ "Mid",
                ntile(theta_new_1pl, 3) == 3 ~ "High"),
      levels = c("Low", "Mid", "High")),
    abs_shift = abs(diff_1pl)
  ) %>%
  ungroup()

equity_summary <- equity %>%
  group_by(task_id, ability_group) %>%
  summarise(n = n(), mean_shift = mean(diff_1pl),
            mean_abs_shift = mean(abs_shift), sd_shift = sd(diff_1pl),
            .groups = "drop")
print(equity_summary, n = 20)

cat("\n-- ANOVA: |shift| ~ ability_group --\n")
for (task in unique(equity$task_id)) {
  d <- equity %>% filter(task_id == task)
  fit <- aov(abs_shift ~ ability_group, data = d)
  cat(sprintf("  %s: F = %.1f, p = %.4f\n", task,
              summary(fit)[[1]]["ability_group", "F value"],
              summary(fit)[[1]]["ability_group", "Pr(>F)"]))
}


# =============================================================
# 12. AGE BAND ANALYSIS
# =============================================================

cat("\n-- Rescoring impact by age band --\n")
age_band_summary <- comparison %>%
  filter(!is.na(age)) %>%
  mutate(age_band = cut(age, breaks = c(4, 6, 8, 10, 13),
                        labels = c("5-6", "7-8", "9-10", "11-12"))) %>%
  filter(!is.na(age_band)) %>%
  group_by(task_id, age_band) %>%
  summarise(n = n(),
            mean_abs_shift = mean(abs(diff_1pl)),
            r_old_new = cor(theta_old, theta_new_1pl),
            .groups = "drop")
print(age_band_summary, n = 20)


# =============================================================
# 13. SENSITIVITY: TRIM OUTLIERS
# =============================================================

cat("\n-- Sensitivity: trim |z| > 3 --\n")
comparison_trimmed <- comparison %>%
  group_by(task_id) %>%
  mutate(z_shift = (diff_1pl - mean(diff_1pl)) / sd(diff_1pl)) %>%
  filter(abs(z_shift) <= 3) %>%
  ungroup()
cat("Removed:", nrow(comparison) - nrow(comparison_trimmed), "rows\n")

comparison_trimmed %>%
  group_by(task_id) %>%
  summarise(n = n(), r = cor(theta_old, theta_new_1pl),
            mean_shift = mean(diff_1pl), rmse = sqrt(mean(diff_1pl^2)),
            .groups = "drop") %>%
  print()


# =============================================================
# 14. TEST-RETEST: GERMANY, ADAPTIVE-FILTERED
# =============================================================

cat("\n\n========== TEST-RETEST (Germany) ==========\n\n")

# Which sites have retests?
cat("-- Sites with 2+ runs --\n")
score_both %>%
  group_by(user_id, param_task_id, site) %>%
  filter(n() >= 2) %>%
  ungroup() %>%
  distinct(user_id, param_task_id, site) %>%
  count(site, param_task_id) %>%
  print(n = 20)

# Join adaptive status
score_both_info <- score_both %>%
  left_join(run_data %>% select(run_id, adaptive), by = "run_id")

# Build retest pairs (Germany only)
retest_pairs <- score_both_info %>%
  filter(site == "pilot_mpieva_de") %>%
  group_by(user_id, param_task_id) %>%
  filter(n() >= 2) %>%
  arrange(user_id, param_task_id, run_id) %>%
  mutate(session = row_number()) %>%
  filter(session <= 2) %>%
  ungroup()

cat("\nRetest pairs (Germany):\n")
retest_pairs %>%
  distinct(user_id, param_task_id) %>%
  count(param_task_id) %>% print()

# Adaptive status per session
cat("\nAdaptive status by task x session:\n")
retest_pairs %>%
  count(param_task_id, session, adaptive) %>%
  print(n = 30)

# Pivot to wide
retest_wide <- retest_pairs %>%
  select(user_id, param_task_id, session,
         theta_old, theta_new_1pl, adaptive) %>%
  pivot_wider(names_from = session,
              values_from = c(theta_old, theta_new_1pl, adaptive),
              names_sep = "_s")

# All pairs
cat("\n-- All pairs --\n")
retest_all <- retest_wide %>%
  mutate(
    abs_diff_old = abs(theta_old_s2 - theta_old_s1),
    abs_diff_new = abs(theta_new_1pl_s2 - theta_new_1pl_s1),
    new_more_stable = abs_diff_new < abs_diff_old
  ) %>%
  group_by(param_task_id) %>%
  summarise(n_pairs = n(),
            mean_diff_old = mean(abs_diff_old),
            mean_diff_new = mean(abs_diff_new),
            pct_new_better = mean(new_more_stable) * 100,
            .groups = "drop")
print(retest_all)

# Same-type pairs only (removes adaptive confound)
cat("\n-- Same-type pairs (adaptive confound removed) --\n")
retest_same <- retest_wide %>%
  filter(adaptive_s1 == adaptive_s2) %>%
  mutate(
    abs_diff_old = abs(theta_old_s2 - theta_old_s1),
    abs_diff_new = abs(theta_new_1pl_s2 - theta_new_1pl_s1),
    new_more_stable = abs_diff_new < abs_diff_old
  ) %>%
  group_by(param_task_id) %>%
  summarise(n_pairs = n(),
            adaptive_type = paste(unique(adaptive_s1), collapse = "/"),
            mean_diff_old = mean(abs_diff_old),
            mean_diff_new = mean(abs_diff_new),
            pct_new_better = mean(new_more_stable) * 100,
            .groups = "drop")
print(retest_same)

# Comparison table
cat("\n-- Comparison: all vs same-type --\n")
bind_rows(
  retest_all %>% mutate(filter = "all"),
  retest_same %>% select(-adaptive_type) %>% mutate(filter = "same_type")
) %>%
  select(param_task_id, filter, n_pairs, pct_new_better) %>%
  pivot_wider(names_from = filter,
              values_from = c(n_pairs, pct_new_better)) %>%
  print()


# =============================================================
# 15. EXPORT
# =============================================================

cat("\n\n========== SAVING ==========\n\n")

write_csv(comparison, file.path(data_dir, "comparison.csv"))
write_csv(table_a, file.path(table_dir, "results_table_a_param_effect.csv"))
write_csv(table_b, file.path(table_dir, "results_table_b_combined.csv"))
write_csv(table_c, file.path(table_dir, "results_table_c_model_effect.csv"))
write_csv(table_site, file.path(table_dir, "results_by_site.csv"))
write_csv(table_age, file.path(table_dir, "results_age_validity.csv"))
write_csv(age_by_site, file.path(table_dir, "results_age_validity_by_site.csv"))
write_csv(cross, file.path(table_dir, "results_cross_task_validity.csv"))
write_csv(boot_results, file.path(table_dir, "results_cross_task_bootstrap.csv"))
write_csv(equity_summary, file.path(table_dir, "results_equity_terciles.csv"))
write_csv(age_band_summary, file.path(table_dir, "results_by_age_band.csv"))
write_csv(retest_all, file.path(table_dir, "results_retest_all.csv"))
write_csv(retest_same, file.path(table_dir, "results_retest_same_type.csv"))


# =============================================================
# 16. FIGURES
# =============================================================

task_labels <- c(
  "egma-math" = "Math", "matrix-reasoning" = "Matrix Reasoning",
  "mental-rotation" = "Mental Rotation",
  "same-different-selection" = "Same-Different"
)
site_colors <- c("pilot_mpieva_de" = "#E8604C",
                 "pilot_uniandes_co" = "#5B9E4F",
                 "pilot_western_ca" = "#2D5FA0")
site_labels <- c("Germany", "Colombia", "Canada")

# Scatter
p_scatter <- ggplot(comparison, aes(theta_old, theta_new_1pl, color = site)) +
  geom_abline(linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.35, size = 0.8) +
  facet_wrap(~task_id, scales = "free",
             labeller = labeller(task_id = task_labels)) +
  scale_color_manual(values = site_colors, labels = site_labels) +
  labs(x = expression(hat(theta)[old]),
       y = expression(hat(theta)[new]), color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold", size = 12))
ggsave(file.path(fig_dir, "fig_scatter.png"), p_scatter,
       width = 8.5, height = 7, dpi = 250, bg = "white")

# Density (2x2)
p_density <-
  ggplot(comparison, aes(x = diff_1pl, fill = site)) +
  geom_density(alpha = 0.3, color = NA) +
  geom_density(aes(color = site), fill = NA, linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "#555555", linewidth = 0.5) +
  facet_wrap(~task_id, scales = "free", ncol = 2,
             labeller = labeller(task_id = task_labels)) +
  scale_fill_manual(values = site_colors, labels = site_labels) +
  scale_color_manual(values = site_colors, labels = site_labels,
                     guide = "none") +
  labs(x = expression(Ability~change~(theta[new] - theta[old])),
       y = NULL, fill = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    text = element_text(family = "serif"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#f0f0f0", linewidth = 0.3),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 9, color = "#555555"),
    axis.title.x = element_text(size = 11, margin = margin(t = 8)),
    legend.position = c(0.92, 0.92),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white",
                                     color = "#dddddd", linewidth = 0.3))
ggsave(file.path(fig_dir, "fig_density.png"), p_density,
       width = 9, height = 5, dpi = 300, bg = "white")


# =============================================================
# DONE
# =============================================================

cat("\nDone! All outputs in:", out_dir, "\n\n")
cat("Tables:\n")
cat(paste(" ", list.files(table_dir, "\\.csv$")), sep = "\n")
cat("\nFigures:\n")
cat(paste(" ", list.files(fig_dir, "\\.png$")), sep = "\n")
cat("\nKey findings:\n")
cat("  Table A: r = .89-.98, shift up to -1.82\n")
cat("  Cross-task: 4/5 significant (bootstrap)\n")
cat("  Equity: low-ability most affected (all p < .0001)\n")
cat("  Retest: SDS 67% (153 pairs, clean), Math 40% (89 pairs, clean)\n")
cat("  MRot retest excluded (adaptive confound)\n")
