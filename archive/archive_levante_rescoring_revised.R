## ============================================================
## LEVANTE CAT Rescoring: REVISED Empirical Analysis
## Fixes: (1) duplicate item_uid, (2) model confound, (3) column clashes
## ============================================================

library(tidyverse)
library(here)

# ── 1. Load data ──────────────────────────────────────────────
setwd("~/Documents/levante-runs/levante-runs")
all_sum_scores <- read_rds(here::here(
  "..", "levante-pilots", "01_fetched_data", "task_data_nested.rds"))

new_scores <- read_csv("all_scores.csv") %>%
  filter(metric_type == "ability",
         task_id %in% c("egma-math", "matrix-reasoning",
                        "mental-rotation", "same-different-selection"))

params <- read_csv("item parameter comparison.csv") %>%
  filter(!is.na(difficulty_old), difficulty_old != "",
         !is.na(difficulty), difficulty != "") %>%
  mutate(b_old = as.numeric(difficulty_old),
         b_new = as.numeric(difficulty),
         param_chance = replace_na(as.numeric(chance_level), 0)) %>%
  select(item_uid, param_task_id = task_id, param_chance, b_old, b_new) %>%
  ## FIX #1: deduplicate params BEFORE joining
  distinct(item_uid, .keep_all = TRUE)

cat("Params after dedup:", nrow(params), "unique items\n")
count(params, param_task_id)

# ── 2. Unnest trials ─────────────────────────────────────────

trials <- all_sum_scores %>% unnest(data)
cat("Trial rows:", nrow(trials), "\n")

# ── 3. EAP scoring function ──────────────────────────────────

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

# ── 4. Merge & score under BOTH old AND new params ───────────
#    FIX #2: Score under new params with same 1PL EAP method
#    This isolates the parameter effect from the model effect

merged <- trials %>%
  filter(!is.na(correct)) %>%
  mutate(correct_num = as.integer(correct)) %>%
  ## FIX #3: params columns already renamed to avoid clashes
  inner_join(params, by = "item_uid") %>%
  filter(param_task_id %in% c("egma-math", "matrix-reasoning",
                               "mental-rotation", "same-different-selection"))

cat("\nMerged trial rows:", nrow(merged), "\n")
cat("Unique children:", n_distinct(merged$user_id), "\n")
cat("By task:\n")
print(count(merged, param_task_id))

# Check: no many-to-many warning should appear now
# If it does, there are still duplicates in trials (multiple runs)
# That's expected and fine -- we group by user_id + run_id below

# Score each child × task under OLD and NEW params (same 1PL EAP method)
score_both <- merged %>%
  group_by(user_id, run_id, param_task_id, site) %>%
  filter(n() >= 5) %>%
  summarise(
    n_items = n(),
    est_old = list(eap_score(correct_num, b_old, c_param = param_chance)),
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

# Keep latest run per child × task (if retests exist)
score_latest <- score_both %>%
  group_by(user_id, param_task_id) %>%
  slice_max(row_number(), n = 1, with_ties = FALSE) %>%
  ungroup()

# ── 5. Also get the multigroup IRT scores (for comparison) ───

score_mg <- new_scores %>%
  group_by(user_id, task_id) %>%
  slice_max(time_started, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(user_id, task_id, age,
         theta_new_mg = metric_value, se_new_mg = metric_se) %>%
  mutate(across(c(theta_new_mg, se_new_mg, age), as.numeric))

# ── 6. Join all three estimates ──────────────────────────────

comparison <- score_latest %>%
  rename(task_id = param_task_id) %>%
  inner_join(score_mg, by = c("user_id", "task_id")) %>%
  mutate(
    diff_1pl = theta_new_1pl - theta_old,   # same method, different params
    diff_mg  = theta_new_mg - theta_old      # different method + different params
  )

cat("\n========== RESULTS ==========\n")
cat("Matched children × tasks:", nrow(comparison), "\n\n")

# ── 7. Summary: isolate parameter effect vs model effect ─────

cat("── A. Parameter effect only (same 1PL EAP, old vs new params) ──\n")
comparison %>%
  group_by(task_id) %>%
  summarise(
    n = n(),
    r_old_new1pl = cor(theta_old, theta_new_1pl),
    mean_shift = mean(diff_1pl),
    sd_shift = sd(diff_1pl),
    rmse = sqrt(mean(diff_1pl^2)),
    .groups = "drop"
  ) %>%
  print()

cat("\n── B. Combined effect (1PL-old vs multigroup-new) ──\n")
comparison %>%
  group_by(task_id) %>%
  summarise(
    n = n(),
    r_old_mg = cor(theta_old, theta_new_mg),
    mean_shift = mean(diff_mg),
    sd_shift = sd(diff_mg),
    rmse = sqrt(mean(diff_mg^2)),
    .groups = "drop"
  ) %>%
  print()

cat("\n── C. Model effect only (same new params, 1PL vs multigroup) ──\n")
comparison %>%
  group_by(task_id) %>%
  summarise(
    n = n(),
    r_1pl_mg = cor(theta_new_1pl, theta_new_mg),
    mean_diff = mean(theta_new_mg - theta_new_1pl),
    sd_diff = sd(theta_new_mg - theta_new_1pl),
    .groups = "drop"
  ) %>%
  print()

# ── 8. By site ───────────────────────────────────────────────

cat("\n── By task × site (parameter effect only) ──\n")
comparison %>%
  group_by(task_id, site) %>%
  summarise(
    n = n(),
    r = cor(theta_old, theta_new_1pl),
    mean_shift = mean(diff_1pl),
    sd_shift = sd(diff_1pl),
    .groups = "drop"
  ) %>%
  print(n = 30)

# ── 9. Validation ────────────────────────────────────────────

cat("\n── Age validity ──\n")
comparison %>%
  filter(!is.na(age)) %>%
  group_by(task_id) %>%
  summarise(
    n = n(),
    r_age_old    = cor(theta_old, age, use = "complete"),
    r_age_new1pl = cor(theta_new_1pl, age, use = "complete"),
    r_age_mg     = cor(theta_new_mg, age, use = "complete"),
    .groups = "drop"
  ) %>%
  print()

cat("\n── Cross-task convergent validity ──\n")
# Pivot to wide
wide <- comparison %>%
  select(user_id, task_id, theta_old, theta_new_1pl, theta_new_mg) %>%
  pivot_wider(names_from = task_id,
              values_from = c(theta_old, theta_new_1pl, theta_new_mg))

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
  m1 <- paste0("theta_new_mg_", pair[1])
  m2 <- paste0("theta_new_mg_", pair[2])

  d <- wide %>% filter(!is.na(.data[[o1]]), !is.na(.data[[o2]]),
                        !is.na(.data[[n1]]), !is.na(.data[[n2]]))
  if (nrow(d) < 10) return(NULL)

  tibble(
    pair = paste0(pair[1], " × ", pair[2]),
    n = nrow(d),
    r_old    = cor(d[[o1]], d[[o2]]),
    r_new1pl = cor(d[[n1]], d[[n2]]),
    r_mg     = cor(d[[m1]], d[[m2]]),
    delta_1pl = cor(d[[n1]], d[[n2]]) - cor(d[[o1]], d[[o2]]),
    delta_mg  = cor(d[[m1]], d[[m2]]) - cor(d[[o1]], d[[o2]])
  )
})
print(cross)

# ── 10. Plots (for poster) ───────────────────────────────────

library(patchwork)

# Use parameter-effect-only comparison for poster (cleaner story)
p1 <- ggplot(comparison, aes(theta_old, theta_new_1pl, color = site)) +
  geom_abline(linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.35, size = 0.8) +
  facet_wrap(~task_id, scales = "free") +
  scale_color_manual(values = c(
    "pilot_mpieva_de" = "#E8604C",
    "pilot_uniandes_co" = "#5B9E4F",
    "pilot_western_ca" = "#2D5FA0"),
    labels = c("Germany", "Colombia", "Canada")) +
  labs(title = "Operational (Old) vs Rescored (New) Ability",
       subtitle = "Same 1PL EAP method, different item parameters",
       x = expression(hat(theta)[old]),
       y = expression(hat(theta)[new]),
       color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold", size = 12))

ggsave("poster_empirical_scatter.png", p1, width = 8.5, height = 7, dpi = 250)
cat("\nSaved poster_empirical_scatter.png\n")

# Score shift density by site
p2 <- ggplot(comparison, aes(diff_1pl, fill = site)) +
  geom_density(alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~task_id, scales = "free") +
  scale_fill_manual(values = c(
    "pilot_mpieva_de" = "#E8604C",
    "pilot_uniandes_co" = "#5B9E4F",
    "pilot_western_ca" = "#2D5FA0"),
    labels = c("Germany", "Colombia", "Canada")) +
  labs(title = "Score Shift by Site (New − Old)",
       subtitle = "Same 1PL EAP method, different item parameters",
       x = "Ability change", fill = NULL) +
  theme_minimal()

ggsave("poster_empirical_density.png", p2, width = 10, height = 6, dpi = 250)
cat("Saved poster_empirical_density.png\n")

cat("\nDone! Check the three summary tables above.\n")
cat("Key question: does Table A (parameter-only) look similar to Table B (combined)?\n")
cat("If yes, the model confound is small and your original numbers hold.\n")
cat("If no, the model difference matters and you should report Table A for the poster.\n")
# Poster-ready density plot: Score Shift by Site
# Run this with your `comparison` dataframe to generate the real figure

library(ggplot2)

# Clean task labels
task_labels <- c(
  "egma-math"                = "Math",
  "matrix-reasoning"         = "Matrix Reasoning",
  "mental-rotation"          = "Mental Rotation",
  "same-different-selection"  = "Same–Different"
)

ggplot(comparison, aes(x = diff_1pl, fill = site)) +
  geom_density(alpha = 0.35, color = NA) +
  geom_density(aes(color = site), fill = NA, linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#555555", linewidth = 0.4) +
  facet_wrap(~task_id, scales = "free", labeller = labeller(task_id = task_labels),ncol = 4) +
  scale_fill_manual(
    values = c(
      "pilot_mpieva_de"   = "#E8604C",
      "pilot_uniandes_co" = "#5B9E4F",
      "pilot_western_ca"  = "#2D5FA0"
    ),
    labels = c("Germany", "Colombia", "Canada")
  ) +
  scale_color_manual(
    values = c(
      "pilot_mpieva_de"   = "#E8604C",
      "pilot_uniandes_co" = "#5B9E4F",
      "pilot_western_ca"  = "#2D5FA0"
    ),
    labels = c("Germany", "Colombia", "Canada"),
    guide = "none"
  ) +
  labs(
    x = expression(Ability~change~(theta[new] - theta[old])),
    y = NULL,
    fill = NULL
  ) +
  theme_minimal(base_size = 10) 

ggsave("poster_density.png", p, width = 7.5, height = 4.5, dpi = 300, bg = "white")
