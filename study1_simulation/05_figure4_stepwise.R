## ============================================================
## Re-knit Figure 4 (stepwise RMSE) with new high-rep data
##
## Inputs:
##   output/data/res_stepwise_v2.rds  (200 reps x 300 examinees,
##                                      from rerun_stepwise_v2.R)
##   output/data/res.rds              (used only to compute Full-100
##                                     reference line at rho = 0.75)
##
## Output:
##   figures/paper/fig_sim_stepwise.png  (overwrites the original)
##   output/tables/results_sim_stepwise.csv
##
## Styling matches the original Figure 4 in run_simulation_analysis_v2.R
## (lines 192-212), with corrected legend mapping so each label has
## a matching colored line swatch.
##
## Usage: run from the repo root (open rescoring-cat.Rproj).
## Runtime: < 5 seconds.
## ============================================================
library(tidyverse)
library(here)   # paths resolve relative to repo root (.here marker)

# -- 1. Read data ----------------------------------------------
stepwise <- readRDS(here("output", "data", "res_stepwise_v2.rds"))
res_raw  <- readRDS(here("output", "data", "res.rds"))
cat("Loaded res_stepwise_v2.rds:", nrow(stepwise), "rows\n")
cat("Conditions:", paste(unique(stepwise$condition), collapse = ", "), "\n")

# -- 2. Compute Full-100 reference at rho = 0.75 ---------------
# Match the original analysis script's computation: mean cell-RMSE for
# Full-100 across all rho cells (Full-100 doesn't depend on rho since
# it scores all 100 items with the true bank).
full100_rmse <- res_raw %>%
  filter(condition == "Full-100", mode == "correlated") %>%
  group_by(rho) %>%
  summarise(cell_rmse = sqrt(mean(mse)), .groups = "drop") %>%
  pull(cell_rmse) %>%
  mean()

cat(sprintf("Full-100 reference RMSE: %.3f\n", full100_rmse))

# -- 3. Verify Figure 4 expectations ---------------------------
cat("\n=== CAT-O step 1-5 (was the bump there before?) ===\n")
stepwise %>%
  filter(condition == "CAT-O", step <= 5) %>%
  select(step, rmse) %>%
  print()

# -- 4. Plot ---------------------------------------------------
# Add Full-100 as a synthetic row so it shows up in the legend with
# a matching swatch, instead of being an unlabeled hline + annotate().
full100_row <- tibble(
  step      = range(stepwise$step),
  rmse      = full100_rmse,
  condition = "Full-100"
)

plot_data <- bind_rows(
  stepwise %>% filter(condition %in% c("CAT-O", "CAT-R", "CAT-B")),
  full100_row
) %>%
  mutate(condition = factor(condition,
                            levels = c("CAT-O", "CAT-R", "CAT-B", "Full-100")))

step_colors <- c(
  "CAT-O"    = "#E8604C",
  "CAT-R"    = "#D4A843",
  "CAT-B"    = "#2D5FA0",
  "Full-100" = "#C48FB5"
)

step_linetypes <- c(
  "CAT-O"    = "solid",
  "CAT-R"    = "solid",
  "CAT-B"    = "solid",
  "Full-100" = "dotted"
)

base_theme <- theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.background = element_rect(fill = "white", color = NA))

p_stepwise <- ggplot(plot_data,
                     aes(x = step, y = rmse,
                         color = condition, linetype = condition)) +
  geom_line(linewidth = 1) +
  # Points only on the simulated conditions, not on the Full-100 ref line.
  geom_point(data = filter(plot_data, condition != "Full-100"),
             size = 1.5) +
  scale_color_manual(values = step_colors,
                     breaks = c("CAT-O", "CAT-R", "CAT-B", "Full-100")) +
  scale_linetype_manual(values = step_linetypes,
                        breaks = c("CAT-O", "CAT-R", "CAT-B", "Full-100")) +
  labs(x = "Number of items administered",
       y = "RMSE",
       color = NULL, linetype = NULL,
       caption = expression("Stepwise RMSE at " * rho * " = 0.75, n = 300, 200 replications")) +
  base_theme

# -- 5. Save ---------------------------------------------------
fig_dir   <- here("figures", "paper")
table_dir <- here("output", "tables")
dir.create(fig_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(fig_dir, "fig_sim_stepwise.png"), p_stepwise,
       width = 7, height = 4.5, dpi = 300, bg = "white")

write_csv(stepwise, file.path(table_dir, "results_sim_stepwise.csv"))

cat("\nWrote:\n")
cat("  ", file.path(fig_dir, "fig_sim_stepwise.png"), "\n")
cat("  ", file.path(table_dir, "results_sim_stepwise.csv"), "\n")

cat("\nNote: paper figure caption should be updated from\n")
cat("  '50 replications' -> '200 replications'\n")
