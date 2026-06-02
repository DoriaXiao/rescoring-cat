## Install all R packages used in this repository.
## Run once: Rscript DEPENDENCIES.R   (or source() it in R)

pkgs <- c(
  "tidyverse",  # dplyr, tidyr, purrr, ggplot2, readr, ...
  "here",       # repo-root-relative paths
  "patchwork",  # multi-panel figures
  "ggrepel",    # non-overlapping labels
  "withr",      # reproducible seeds in stepwise plot
  "tibble",
  "mirt"        # item calibration (study1/03 N->rho; multigroup IRT)
)

to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install)) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  message("All required packages already installed.")
}
