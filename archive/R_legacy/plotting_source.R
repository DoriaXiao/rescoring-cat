plot_metric_accuracy_v4 <- function(
    results_df,
    metric = "rmse",
    methods = NULL,            # optional subset of labeled methods to display
    family  = c("auto","correlated","additive","all"),
    xvar    = c("auto","rho","sigma"),
    error_bar   = TRUE,
    show_labels = TRUE,
    # smoothing
    smooth = FALSE,
    smooth_method = c("loess","gam","none"),
    smooth_span   = 0.5,
    smooth_degree = 2,
    smooth_se     = FALSE,
    show_raw_lines = TRUE
){
  library(dplyr); library(ggplot2); library(ggrepel); library(tidyr)
  
  family <- match.arg(family)
  smooth_method <- match.arg(smooth_method)
  xvar <- match.arg(xvar)
  
  # --- palettes
  okabe_ito <- c(
    black="#000000", orange="#E69F00", sky="#56B4E9", bluish="#009E73",
    yellow="#F0E442", blue="#0072B2", vermillion="#D55E00", purple="#CC79A7"
  )
  
  # --- label mapping for lines (matches your meeting terms)
  label_design <- function(design, select_bank, score_with) {
    dplyr::case_when(
      design == "CAT"   & select_bank == "bhat_cat"  & score_with == "noisy" ~ "CAT (operational)",
      design == "CAT"   & select_bank == "bhat_cat"  & score_with == "true"  ~ "CAT (rescoring)",
      design == "CAT"   & select_bank == "bhat_true" & score_with == "true"  ~ "CAT (benchmark)",
      design == "CAT"   & select_bank == "b"                                 ~ "CAT (oracle)",
      design == "Fixed" & score_with == "true" ~ "Fixed (20)",
      design == "Full"  & score_with == "true" ~ "Full (100)",
      TRUE ~ NA_character_
    )
  }
  
  # --- choose family + xvar
  df <- results_df
  if (family == "correlated") df <- df %>% filter(mode == "correlated")
  if (family == "additive")   df <- df %>% filter(mode == "additive")
  if (xvar == "auto") {
    xvar <- if (n_distinct(df$mode) == 1 && unique(df$mode) == "correlated") "rho" else "sigma"
  }
  
  # --- bias condition labels (including correlated family)
  BiasCond <- function(mode, bias_type, b_bias) {
    dplyr::case_when(
      mode == "correlated"                         ~ "Correlated (ρ)",
      bias_type == "random"                        ~ "Random Error",
      bias_type == "structured" & b_bias <  0     ~ "Structured: Underestimation (-0.3)",
      bias_type == "structured" & b_bias ==  0.3     ~ "Structured: Overestimation (+0.3)",
      bias_type == "structured" & b_bias ==1.0     ~ "Structured: Overestimation (+1)",
      TRUE                                         ~ "Other"
    )
  }
  
  metric <- tolower(metric)
  df <- df %>%
    mutate(
      Design = label_design(design, select_bank, score_with),
      BiasCond = BiasCond(mode, bias_type, b_bias),
      Shift = ifelse(is.na(b_bias) | b_bias == 0, "Δb=0", paste0("Δb=", b_bias)),
      xval = if (xvar == "rho") rho else sigma
    ) %>%
    filter(!is.na(Design), !is.na(xval))
  
  if (!is.null(methods)) df <- df %>% filter(Design %in% methods)
  
  # --- aggregate over replicates
  plot_data <- df %>%
    group_by(BiasCond, Shift, Design, xval) %>%
    summarise(
      metric_val = mean(.data[[metric]], na.rm = TRUE),
      metric_se  = sd(.data[[metric]], na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    filter(!is.na(metric_val))
  
  if (!nrow(plot_data)) stop("No rows to plot after filtering/aggregation.")
  
  # --- color map (Okabe-Ito, deterministic)
  base_cols <- c(
    "CAT (operational)" = okabe_ito[["vermillion"]],
    "CAT (rescoring)"   = okabe_ito[["orange"]],
    "CAT (benchmark)"   = okabe_ito[["blue"]],
    "CAT (oracle)"      = okabe_ito[["sky"]],
    "Fixed (20)"        = okabe_ito[["black"]],
    "Full (100)"        = okabe_ito[["purple"]]
  )
  present <- unique(plot_data$Design)
  col_map <- base_cols[names(base_cols) %in% present]
  if (length(col_map) < length(present)) {
    extra <- setdiff(present, names(col_map))
    pool  <- rep(unname(okabe_ito), length.out = length(extra))
    names(pool) <- extra
    col_map <- c(col_map, pool)
  }
  
  # --- axis/labels
  y_label <- dplyr::case_when(
    metric == "rmse"               ~ "Root Mean Square Error (RMSE)",
    metric == "bias"               ~ "Average Estimation Bias",
    metric == "pearson_r"          ~ "Pearson Correlation (θ̂, θ)",
    metric == "kendall_tau"        ~ "Kendall’s Tau (θ̂, θ)",
    metric == "spearman_rho"       ~ "Spearman’s Rho (θ̂, θ)",
    metric == "targeting_distance" ~ "Targeting Distance (mean |θ − b|)",
    metric == "targeting_distance_b" ~ "Targeting: TRUE bank mean |θ − b|",
    metric == "targeting_distance_bprime" ~ "Targeting: admin bank mean |θ − b′|",
    TRUE ~ metric
  )
  x_label <- if (xvar == "rho") expression(rho) else expression(sigma)
  title_label <- if (xvar == "rho") "Performance vs Alignment ρ" else "Performance vs Calibration Noise σ"
  subtitle_label <- "Error bars = ±1 SE"
  
  # --- label positions
  lab_df <- plot_data %>%
    group_by(BiasCond, Shift, Design) %>%
    slice_max(xval, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # --- build plot
  g <- ggplot(plot_data, aes(x = xval, y = metric_val, color = Design, group = Design))
  
  if (smooth && smooth_method == "loess") {
    if (show_raw_lines) g <- g + geom_line(linewidth = 0.6, alpha = 0.35)
    g <- g + geom_smooth(
      method = "loess", se = smooth_se, span = smooth_span,
      method.args = list(degree = smooth_degree, family = "gaussian"),
      linewidth = 1.05
    )
  } else if (smooth && smooth_method == "gam") {
    if (show_raw_lines) g <- g + geom_line(linewidth = 0.6, alpha = 0.35)
    g <- g + geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = smooth_se, linewidth = 1.05)
  } else {
    g <- g + geom_line(linewidth = 1.05, alpha = 0.95) + geom_point(size = 2.3)
  }
  
  if (isTRUE(error_bar)) {
    g <- g + geom_errorbar(aes(ymin = metric_val - metric_se, ymax = metric_val + metric_se),
                           width = 0.02, alpha = 0.7)
  }
  if (metric == "bias") g <- g + geom_hline(yintercept = 0, linetype = "dashed", color = "gray60")
  
  if (isTRUE(show_labels) && nrow(lab_df)) {
    g <- g + ggrepel::geom_text_repel(
      data = lab_df, aes(label = Design), size = 3.4,
      segment.alpha = 0.25, direction = "y", box.padding = 0.3, max.overlaps = Inf
    )
  }
  
  # facet by bias condition; for additive this splits Random vs Structured ±
  g +
    facet_grid(cols = vars(BiasCond)) +
    scale_color_manual(values = col_map, breaks = present) +
    labs(title = title_label, subtitle = subtitle_label, x = x_label, y = y_label, color = NULL) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold", size = 12))
}
