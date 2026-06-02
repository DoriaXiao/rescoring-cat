# Hardened helper: replace your existing add_calibration_error() with this
add_calibration_error <- function(
    items_core,
    sigma,
    bias_type       = "random",
    b_bias          = 0,
    err_mode        = c("additive","correlated"),  # was `mode`
    rho             = 1.0,
    standardize     = FALSE
){
  err_mode <- match.arg(err_mode, choices = c("additive","correlated"))
  bias_type <- as.character(bias_type)[1]
  
  items <- items_core
  
  # Always keep originals
  a_true <- items$a_true
  b_true <- items$b_true
  
  if (err_mode == "additive") {
    # global shift + random noise
    b_prime <- b_true + b_bias + rnorm(length(b_true), mean = 0, sd = sigma)
  } else {
    # correlated cross-population misalignment:
    # b' = rho * b + sqrt(1 - rho^2) * z, z ~ N(0, var(b_true))
    sd_b <- stats::sd(b_true)
    z    <- rnorm(length(b_true), 0, sd_b)
    b_prime <- rho * b_true + sqrt(pmax(0, 1 - rho^2)) * z + b_bias
    # optional random noise on top
    if (sigma > 0) b_prime <- b_prime + rnorm(length(b_prime), 0, sigma)
  }
  
  if (isTRUE(standardize)) {
    # Center/scale b' to match b distribution (if you’ve been doing this)
    b_prime <- (b_prime - mean(b_prime)) / sd(b_prime) * sd(b_true) + mean(b_true)
  }
  
  # If you estimate a' (discrimination) separately under noisy bank, do it here;
  # otherwise keep a_cat = a_true and b_cat = b_prime as the calibration bank.
  items$a_cat <- a_true
  items$b_cat <- b_prime
  
  # If you also keep ML/EAP-estimated versions (hat banks), leave NULL here;
  # your estimate_banks() will fill `bh_cat`, `bh_true`, etc.
  items
}


# --- Rep-aware, field-selectable summarizer -----------------------------------
summarize_rmse_by_step <- function(tr, rep_id = NA_integer_, field = "theta_hat_step_true") {
  est_mat <- tr[[field]]
  steps   <- ncol(est_mat)
  tibble::tibble(
    step       = rep(seq_len(steps), each = length(tr$theta_true)),
    theta_true = rep(tr$theta_true,    times = steps),
    est        = as.vector(est_mat),
    rep        = rep_id
  ) |>
    dplyr::group_by(step, rep) |>
    dplyr::summarise(
      rmse = sqrt(mean((est - theta_true)^2, na.rm = TRUE)),
      .groups = "drop"
    )
}

# --- CAT trace with deterministic max-info selection --------------------------
run_cat_trace_select <- function(
    sigma, b_sd, b_bias = 0, n_persons = 400, max_steps = 30,
    bias_type = "random", sem_threshold = 0.4, min_steps = 5,
    top_k = 1,
    select_bank = c("bhat_cat","bhat_true","bprime","b"),
    score_with  = c("noisy","true"),
    items_pool  = NULL,
    err_mode    = c("additive","correlated"),
    rho  = 1.0, standardize_bprime = FALSE,
    N_cat = 150, N_rescore = 300,
    theta_true = NULL,
    do_stop    = FALSE
){
  select_bank <- match.arg(select_bank, choices = c("bhat_cat","bhat_true","bprime","b"))
  score_with  <- match.arg(score_with,  choices = c("noisy","true"))
  err_mode    <- match.arg(err_mode,    choices = c("additive","correlated"))
  
  items_core <- if (is.null(items_pool)) generate_items_pool(100, b_sd = b_sd) else items_pool
  
  # err_mode passed as the 'mode' argument your helper expects
  items <- add_calibration_error(items_core, sigma, bias_type, b_bias, err_mode, rho, standardize_bprime)
  items <- estimate_banks(items, N_cat = N_cat, N_rescore = N_rescore)
  
  a_sel <- switch(select_bank,
                  "bhat_cat"  = items$a_cat,
                  "bhat_true" = items$a_true,
                  "bprime"    = items$a_cat,
                  "b"         = items$a_true)
  b_sel <- switch(select_bank,
                  "bhat_cat"  = items$bh_cat,
                  "bhat_true" = items$bh_true,
                  "bprime"    = items$b_cat,
                  "b"         = items$b_true)
  
  a_scr <- if (score_with == "noisy") a_sel else items$a_true
  b_scr <- if (score_with == "noisy") b_sel else items$b_true
  
  if (is.null(theta_true)) theta_true <- rnorm(n_persons)
  
  theta_hat_step       <- matrix(NA_real_, nrow = n_persons, ncol = max_steps)
  theta_hat_step_sel   <- matrix(NA_real_, nrow = n_persons, ncol = max_steps)
  theta_hat_step_true  <- matrix(NA_real_, nrow = n_persons, ncol = max_steps)
  
  for (i in 1:n_persons) {
    administered <- integer(0); responses <- integer(0); theta_hat <- 0
    for (step in 1:max_steps) {
      avail <- setdiff(seq_len(nrow(items)), administered)
      if (!length(avail)) break
      
      p   <- plogis(a_sel[avail] * (theta_hat - b_sel[avail]))
      inf <- a_sel[avail]^2 * p * (1 - p)
      
      pick      <- avail[order(inf, decreasing = TRUE)[1:min(top_k, length(inf))]]
      next_item <- pick[1]  # deterministic
      
      p_true     <- plogis(items$a_true[next_item] * (theta_true[i] - items$b_true[next_item]))
      responses  <- c(responses, rbinom(1, 1, p_true))
      administered <- c(administered, next_item)
      
      theta_hat <- estimate_theta(a_sel[administered], b_sel[administered], responses)
      theta_hat_step_sel [i, step] <- theta_hat
      theta_hat_step_true[i, step] <- estimate_theta(items$a_true[administered], items$b_true[administered], responses)
      theta_hat_step     [i, step] <- estimate_theta(a_scr[administered],        b_scr[administered],        responses)
      
      if (do_stop && step >= min_steps) {
        p_now <- plogis(a_sel[administered] * (theta_hat - b_sel[administered]))
        info  <- sum(a_sel[administered]^2 * p_now * (1 - p_now))
        if (info > 0 && 1/sqrt(info) < sem_threshold) break
      }
    }
  }
  
  list(
    theta_true          = theta_true,
    theta_hat_step      = theta_hat_step,       # scored with (a_scr, b_scr)
    theta_hat_step_sel  = theta_hat_step_sel,   # selection-only estimate
    theta_hat_step_true = theta_hat_step_true,  # rescored on TRUE bank
    max_steps           = max_steps,
    meta = list(bias_type=bias_type, err_mode=err_mode, rho=rho, sigma=sigma, b_bias=b_bias)
  )
}

# --- Full-form trace: stepwise scoring under b (TRUE) -------------------------
run_full_form_trace <- function(
    sigma, b_sd, b_bias = 0,
    n_persons = 400, n_items = 100,
    bias_type = "random",
    items_pool = NULL,
    err_mode   = c("additive","correlated"),
    rho = 1.0, standardize_bprime = FALSE,
    theta_true = NULL
){
  err_mode <- match.arg(err_mode, choices = c("additive","correlated"))
  items_core <- if (is.null(items_pool)) generate_items_pool(n_items, b_sd = b_sd) else items_pool
  items <- add_calibration_error(items_core, sigma, bias_type, b_bias, err_mode, rho, standardize_bprime)
  
  if (is.null(theta_true)) theta_true <- rnorm(n_persons)
  P_true <- plogis(outer(theta_true, items$b_true, "-") *
                     matrix(items$a_true, n_persons, n_items, byrow = TRUE))
  R <- matrix(rbinom(n_persons * n_items, 1, as.vector(P_true)), n_persons, n_items)
  
  theta_hat_step_true <- matrix(NA_real_, n_persons, n_items)
  for (i in 1:n_persons) for (k in 1:n_items) {
    use <- 1:k
    theta_hat_step_true[i, k] <- estimate_theta(items$a_true[use], items$b_true[use], R[i, use])
  }
  
  list(
    theta_true = theta_true,
    theta_hat_step_true = theta_hat_step_true,
    max_steps = n_items,
    meta = list(bias_type=bias_type, err_mode=err_mode, rho=rho, sigma=sigma, b_bias=b_bias)
  )
}

# --- Fixed-form trace: stepwise scoring under b (TRUE) ------------------------
run_fixed_trace <- function(
    sigma, b_sd, b_bias = 0,
    test_length = 20, n_persons = 400,
    bias_type = "random",
    items_pool = NULL, fixed_indices = NULL,
    err_mode = c("additive","correlated"),
    rho = 1.0, standardize_bprime = FALSE,
    theta_true = NULL
){
  err_mode <- match.arg(err_mode, choices = c("additive","correlated"))
  items_core <- if (is.null(items_pool)) generate_items_pool(100, b_sd = b_sd) else items_pool
  items <- add_calibration_error(items_core, sigma, bias_type, b_bias, err_mode, rho, standardize_bprime)
  
  if (is.null(fixed_indices)) fixed_indices <- order(-items_core$a_true)[1:test_length]
  if (is.null(theta_true)) theta_true <- rnorm(n_persons)
  
  theta_hat_step_true <- matrix(NA_real_, n_persons, test_length)
  for (i in 1:n_persons) {
    p_true <- plogis(items$a_true[fixed_indices] * (theta_true[i] - items$b_true[fixed_indices]))
    resp   <- rbinom(test_length, 1, p_true)
    for (k in 1:test_length) {
      use <- fixed_indices[1:k]
      theta_hat_step_true[i, k] <- estimate_theta(items$a_true[use], items$b_true[use], resp[1:k])
    }
  }
  
  list(
    theta_true = theta_true,
    theta_hat_step_true = theta_hat_step_true,
    max_steps = test_length,
    design = "Fixed",
    meta = list(bias_type=bias_type, err_mode=err_mode, rho=rho, sigma=sigma, b_bias=b_bias,
                fixed_indices=fixed_indices)
  )
}

plot_traj_simple <- function(
    res,
    bias_type    = "none",
    b_bias       = 0,
    n_persons    = 500,
    max_steps    = 25,
    fixed_length = 20,
    methods      = c("cat_ops","cat_rescore","full","fixed_point","fixed_traj"),
    err_mode     = c("additive","correlated"),
    rho          = NA_real_,
    standardize_bprime = FALSE,
    items_pool   = NULL,
    fixed_indices= NULL,
    N_cat        = 150,
    N_rescore    = 300,
    seed         = 20250831,
    n_reps       = 40,
    show_ribbon  = TRUE,
    keep_only    = c("cat_ops","cat_rescore","fixed_traj","full")
){
  err_mode <- match.arg(err_mode, choices = c("additive","correlated"))
  
  suppressPackageStartupMessages({
    library(dplyr); library(ggplot2); library(tibble); library(purrr)
  })
  
  sub <- res %>% dplyr::filter(.data$bias_type == bias_type, .data$b_bias == b_bias)
  sigma_use <- if (nrow(sub)) stats::median(sub$sigma, na.rm = TRUE) else 0
  b_sd_use  <- if (nrow(sub)) sub$b_sd[which(!is.na(sub$b_sd))[1]] else 0.75
  if (!is.finite(b_sd_use)) b_sd_use <- 0.75
  rho_use   <- if (is.finite(rho)) rho else if ("rho" %in% names(sub) && any(is.finite(sub$rho))) stats::median(sub$rho, na.rm = TRUE) else 0.85
  
  if (is.null(items_pool))   items_pool   <- generate_items_pool(n_items = 100, b_sd = b_sd_use)
  if (is.null(fixed_indices)) fixed_indices <- order(-items_pool$a_true)[1:fixed_length]
  
  label_map <- c(
    cat_ops     = "CAT (operational)",
    cat_rescore = "CAT (rescoring)",
    full        = "Full-100",
    fixed_point = "Fixed",
    fixed_traj  = sprintf("Fixed-%d", fixed_length)
  )
  
  run_one <- function(rep_idx) {
    withr::with_seed(seed + rep_idx, {
      parts <- list()
      theta_draw <- rnorm(n_persons)
      
      if ("cat_ops" %in% methods) {
        tr <- run_cat_trace_select(
          sigma = sigma_use, b_sd = b_sd_use, b_bias = b_bias,
          n_persons = n_persons, max_steps = max_steps, bias_type = bias_type,
          select_bank = "bhat_cat", score_with = "noisy",
          items_pool = items_pool, err_mode = err_mode, rho = rho_use,
          standardize_bprime = standardize_bprime, N_cat = N_cat, N_rescore = N_rescore,
          theta_true = theta_draw, do_stop = FALSE, top_k = 1
        )
        parts$cat_ops <- summarize_rmse_by_step(tr, rep_id = rep_idx, field = "theta_hat_step") %>%
          dplyr::mutate(method = label_map[["cat_ops"]])
      }
      
      if ("cat_rescore" %in% methods) {
        tr <- run_cat_trace_select(
          sigma = sigma_use, b_sd = b_sd_use, b_bias = b_bias,
          n_persons = n_persons, max_steps = max_steps, bias_type = bias_type,
          select_bank = "bhat_cat", score_with = "true",
          items_pool = items_pool, err_mode = err_mode, rho = rho_use,
          standardize_bprime = standardize_bprime, N_cat = N_cat, N_rescore = N_rescore,
          theta_true = theta_draw, do_stop = FALSE, top_k = 1
        )
        parts$cat_rescore <- summarize_rmse_by_step(tr, rep_id = rep_idx, field = "theta_hat_step_true") %>%
          dplyr::mutate(method = label_map[["cat_rescore"]])
      }
      
      if ("full" %in% methods) {
        tr <- run_full_form_trace(
          sigma = sigma_use, b_sd = b_sd_use, b_bias = b_bias,
          n_persons = n_persons, n_items = 100, bias_type = bias_type,
          items_pool = items_pool, err_mode = err_mode, rho = rho_use,
          standardize_bprime = standardize_bprime, theta_true = theta_draw
        )
        parts$full <- summarize_rmse_by_step(tr, rep_id = rep_idx, field = "theta_hat_step_true") %>%
          dplyr::mutate(method = label_map[["full"]])
      }
      
      if ("fixed_traj" %in% methods) {
        tr <- run_fixed_trace(
          sigma = sigma_use, b_sd = b_sd_use, b_bias = b_bias,
          test_length = fixed_length, n_persons = n_persons, bias_type = bias_type,
          items_pool = items_pool, fixed_indices = fixed_indices,
          err_mode = err_mode, rho = rho_use, standardize_bprime = standardize_bprime,
          theta_true = theta_draw
        )
        parts$fixed_traj <- summarize_rmse_by_step(tr, rep_id = rep_idx, field = "theta_hat_step_true") %>%
          dplyr::mutate(method = label_map[["fixed_traj"]])
      }
      
      dplyr::bind_rows(parts)
    })
  }
  
  df_all <- purrr::map_dfr(seq_len(n_reps), run_one) %>%
    dplyr::filter(.data$method %in% label_map[keep_only])
  
  df_sum <- df_all %>%
    dplyr::group_by(method, step) %>%
    dplyr::summarise(
      mean_rmse = mean(rmse, na.rm = TRUE),
      se        = stats::sd(rmse, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups   = "drop"
    )
  
  endpts <- df_sum %>%
    dplyr::group_by(method) %>%
    dplyr::filter(step == max(step)) %>%
    dplyr::slice_tail(n = 1)
  
  full_row <- df_sum %>% dplyr::filter(method == "Full-100", step == max(step))
  
  title_txt <- if (err_mode == "correlated") {
    sprintf("Length trajectories under correlated misalignment (ρ = %.2f)", rho_use)
  } else {
    sprintf("Length trajectories under structured shift (Δb = %s)", b_bias)
  }
  # 1) grab the fixed endpoint (a single dot)
  fixed_end <- df_sum %>%
    dplyr::filter(method == sprintf("Fixed-%d", fixed_length), step == fixed_length)
  
  # 2) grab the CAT lines only (no fixed/full trajectories)
  df_lines <- df_sum %>%
    dplyr::filter(method %in% c("CAT (operational)", "CAT (rescoring)"))
  # --- title ---
  title_txt <- if (err_mode == "correlated") {
    sprintf("Length trajectories under correlated misalignment (ρ = %.2f)", rho_use)
  } else {
    sprintf("Length trajectories under structured shift (Δb = %s)", b_bias)
  }
  
  # --- data splits for plotting ---
  # CAT curves (lines)
  df_lines <- df_sum %>%
    dplyr::filter(method %in% c("CAT (operational)", "CAT (rescoring)"))
  
  # Full-100 ceiling (final point)
  full_row <- df_sum %>% dplyr::filter(method == "Full-100", step == max(step))
  
  # Fixed trajectory as dots across all steps
  df_fixed_pts <- df_sum %>%
    dplyr::filter(method == sprintf("Fixed-%d", fixed_length)) %>%
    dplyr::arrange(step)
  
  # Fixed endpoint (for label)
  fixed_end <- df_fixed_pts %>% dplyr::filter(step == fixed_length)
  
  # --- toggle for SE whiskers on fixed dots ---
  fixed_dots_se <- TRUE  # set FALSE if you don't want whiskers
  
  # --- plot ---
  p <- ggplot(df_lines, aes(step, mean_rmse, color = method, fill = method)) +
    (if (show_ribbon) geom_ribbon(aes(ymin = mean_rmse - se, ymax = mean_rmse + se),
                                  alpha = 0.12, color = NA) else NULL) +
    geom_line(linewidth = 1) +
    # Full-100 ceiling (horizontal line + label)
    (if (nrow(full_row)) geom_hline(yintercept = full_row$mean_rmse[1], linetype = 3) else NULL) +
    (if (nrow(full_row)) annotate("text", x = max_steps - 1, y = full_row$mean_rmse[1],
                                  label = "Full-100 ceiling", hjust = 1, vjust = -0.6, size = 3) else NULL) +
    # Fixed trajectory as dots at every k
    (if (nrow(df_fixed_pts)) geom_point(data = df_fixed_pts,
                                        aes(step, mean_rmse),
                                        inherit.aes = FALSE,
                                        size = 1.8, shape = 16, alpha = 0.9,
                                        color = "steelblue3") else NULL) +
    # Optional SE whiskers on fixed dots
    (if (nrow(df_fixed_pts) && fixed_dots_se)
      geom_errorbar(data = df_fixed_pts,
                    aes(x = step, ymin = mean_rmse - se, ymax = mean_rmse + se),
                    inherit.aes = FALSE, width = 0.35, alpha = 0.6,
                    color = "steelblue3") else NULL) +
    # Label the fixed endpoint
    (if (nrow(fixed_end)) annotate("text", x = fixed_length, y = fixed_end$mean_rmse[1],
                                   label = sprintf("Fixed-%d", fixed_length),
                                   vjust = -1.0, hjust = 0.5, size = 3, color = "steelblue3") else NULL) +
    scale_x_continuous(breaks = seq(0, max_steps, by = 5), limits = c(0, max_steps)) +
    labs(x = "# items (k)", y = "RMSE",
         title = title_txt,
         subtitle = "CAT lines (with SE bands); Fixed form as dots at each k; Full-100 is the ceiling",
         color = NULL, fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  
  p
}

