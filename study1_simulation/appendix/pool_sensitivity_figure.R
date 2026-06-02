## ============================================================
## Sensitivity Figure for Appendix
##
## Reads:   output/data/res_sensitivity.rds
## Outputs: figures/paper/fig_sim_sensitivity.png
##
## Shows CAT-O and CAT-R mean curves across rho, with 10th-90th
## percentile bands across 10 pools, to demonstrate that the
## rescoring effect is robust to pool draw.
## ============================================================

library(tidyverse)

sens <- readRDS("output/data/res_sensitivity.rds")

fig_dir <- "figures/paper"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Summarise: mean and percentile bands across pools for each (rho, condition)
summary_bands <- sens %>%
  filter(condition %in% c("CAT-O", "CAT-R", "CAT-B")) %>%
  group_by(rho, condition) %>%
  summarise(
    mean_rmse = mean(rmse),
    lo        = quantile(rmse, 0.10),
    hi        = quantile(rmse, 0.90),
    .groups   = "drop"
  )

cond_colors <- c("CAT-O" = "#E8604C", "CAT-R" = "#D4A843",
                 "CAT-B" = "#2D5FA0")

p_sens <- ggplot(summary_bands,
                 aes(x = rho, color = condition, fill = condition)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.20, color = NA) +
  geom_line(aes(y = mean_rmse), linewidth = 1) +
  geom_point(aes(y = mean_rmse), size = 1.5) +
  scale_color_manual(values = cond_colors,
                     breaks = c("CAT-O", "CAT-R", "CAT-B")) +
  scale_fill_manual(values = cond_colors,
                    breaks = c("CAT-O", "CAT-R", "CAT-B"),
                    guide = "none") +
  labs(x = expression("Cross-population alignment ("*rho*")"),
       y = "RMSE",
       color = NULL,
       caption = "Lines: mean across 10 pools. Shaded bands: 10th-90th percentile across pools.") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.background = element_rect(fill = "white", color = NA))

ggsave(file.path(fig_dir, "fig_sim_sensitivity.png"), p_sens,
       width = 7, height = 4.5, dpi = 300, bg = "white")

cat("Saved:", file.path(fig_dir, "fig_sim_sensitivity.png"), "\n")
