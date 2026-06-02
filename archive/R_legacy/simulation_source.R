# --- 1. Setup and Core Functions ---

# Libraries 
library(dplyr); library(tidyr); library(purrr); library(ggplot2)
library(mirt); library(mirtCAT); library(patchwork)
# ---- Helpers ---------------------------------------------------------------
# ---- Bank estimation helpers ----
# If you already have your own `estimate_bank()`, keep it; otherwise this minimal
# fallback adds finite-sample noise around b to mimic calibration.
if (!exists("estimate_bank", mode = "function")) {
  estimate_bank <- function(b_vec, N = 300) {
    se <- stats::sd(b_vec, na.rm = TRUE) / sqrt(max(1, N))
    b_vec + stats::rnorm(length(b_vec), mean = 0, sd = se)
  }
}

estimate_banks <- function(items, N_cat = 150, N_rescore = 300) {
  stopifnot(all(c("b_true","b_cat") %in% names(items)))
  items$bh_true <- estimate_bank(items$b_true, N = N_rescore)  # \hat b
  items$bh_cat  <- estimate_bank(items$b_cat,  N = N_cat)      # \hat b′
  items
}
.draw_param <- function(draw = c("vector","range"),
                        vec = NULL, range = c(0,1), n = 20,
                        clamp = NULL) {
  draw <- match.arg(draw)
  vals <- if (draw == "vector") {
    if (is.null(vec) || !length(vec)) stop("Provide a non-empty vector for 'vec'.")
    vec
  } else {
    runif(n, min = range[1], max = range[2])
  }
  if (!is.null(clamp)) {
    stopifnot(length(clamp) == 2)
    vals <- pmin(pmax(vals, clamp[1]), clamp[2])
  }
  vals
}

.validate_select_bank_opts <- function(select_bank, select_bank_opts) {
  opts <- if (!is.null(select_bank_opts)) select_bank_opts else select_bank
  stopifnot(all(opts %in% c("bhat_cat","bhat_true","bprime","b")))
  opts
}

# ---- Build the full parameter grid (no sims) -------------------------------

build_param_grid <- function(
    # σ controls (now symmetric with ρ)
  sigma_draw  = c("range","vector"),
  sigma_range = c(0, 1),
  n_sigma     = 20,
  sigma_vec   = NULL,
  
  # ρ controls
  rho_draw    = c("vector","range"),
  rho_vec     = NULL,
  rho_range   = c(0, 1),
  n_rho       = 20,
  
  # design axes
  b_sd_vec = c(0.5, 0.75, 1.0),
  stopping_rules = c("SEM","fixed","SEM_topk","full"),
  b_bias_vec = c(-0.5, 0.5),
  modes = c("correlated"),
  select_bank = "bhat_cat",
  select_bank_opts = NULL,
  
  # NEW: ignore σ in correlated rows
  ignore_sigma_in_correlated = TRUE,
  sigma_anchor_for_correlated = 0,
  
  # reproducibility
  grid_seed = NULL
){
  if (!is.null(grid_seed)) set.seed(grid_seed)
  
  sigma_draw <- match.arg(sigma_draw)
  rho_draw   <- match.arg(rho_draw)
  
  # draws
  sigma_vals <- .draw_param(sigma_draw, sigma_vec, sigma_range, n_sigma,
                            clamp = c(0, Inf))
  rho_vals   <- .draw_param(rho_draw,   rho_vec,   rho_range,   n_rho,
                            clamp = c(-0.999, 0.999))
  
  cat_rules   <- setdiff(stopping_rules, c("fixed","full"))
  other_rules <- intersect(stopping_rules, c("fixed","full"))
  sel_vec     <- .validate_select_bank_opts(select_bank, select_bank_opts)
  
  grids <- list()
  
  # --- CORRELATED (bias_type = 'none') ---
  if ("correlated" %in% modes) {
    sig_cor <- if (ignore_sigma_in_correlated) sigma_anchor_for_correlated else sigma_vals
    
    correlated_cat <- expand.grid(
      sigma = sig_cor, b_sd = b_sd_vec, stopping_rule = cat_rules,
      bias_type = "none", b_bias = 0,
      rho = rho_vals, select_bank = sel_vec, mode = "correlated",
      stringsAsFactors = FALSE
    )
    correlated_other <- expand.grid(
      sigma = sig_cor, b_sd = b_sd_vec, stopping_rule = other_rules,
      bias_type = "none", b_bias = 0,
      rho = rho_vals, select_bank = NA_character_, mode = "correlated",
      stringsAsFactors = FALSE
    )
    grids <- c(grids, list(correlated_cat, correlated_other))
  }
  
  # --- ADDITIVE (random + structured) ---
  if ("additive" %in% modes) {
    add_random_cat <- expand.grid(
      sigma = sigma_vals, b_sd = b_sd_vec, stopping_rule = cat_rules,
      bias_type = "random", b_bias = 0,
      rho = NA_real_, select_bank = sel_vec, mode = "additive",
      stringsAsFactors = FALSE
    )
    add_random_other <- expand.grid(
      sigma = sigma_vals, b_sd = b_sd_vec, stopping_rule = other_rules,
      bias_type = "random", b_bias = 0,
      rho = NA_real_, select_bank = NA_character_, mode = "additive",
      stringsAsFactors = FALSE
    )
    add_struct_cat <- expand.grid(
      sigma = sigma_vals, b_sd = b_sd_vec, stopping_rule = cat_rules,
      bias_type = "structured", b_bias = b_bias_vec,
      rho = NA_real_, select_bank = sel_vec, mode = "additive",
      stringsAsFactors = FALSE
    )
    add_struct_other <- expand.grid(
      sigma = sigma_vals, b_sd = b_sd_vec, stopping_rule = other_rules,
      bias_type = "structured", b_bias = b_bias_vec,
      rho = NA_real_, select_bank = NA_character_, mode = "additive",
      stringsAsFactors = FALSE
    )
    grids <- c(grids, list(add_random_cat, add_random_other, add_struct_cat, add_struct_other))
  }
  
  dplyr::bind_rows(grids)
}

# ---- Helpers ---------------------------------------------------------------

.validate_select_bank_opts <- function(select_bank, select_bank_opts) {
  opts <- if (!is.null(select_bank_opts)) select_bank_opts else select_bank
  stopifnot(all(opts %in% c("bhat_cat","bhat_true","bprime","b")))
  opts
}


# --- Safe helpers -------------------------------------------------------------
safe_cor <- function(x, y, method = c("pearson","spearman","kendall")) {
  method <- match.arg(method)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) >= 3) cor(x[ok], y[ok], method = method) else NA_real_
}
# --- NEW: utilities -----------------------------------------------------------

# Make x have the same mean/SD as ref (if we keep the additive b' path)
standardize_like <- function(x, ref) {
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE) * sd(ref, na.rm=TRUE) + mean(ref, na.rm=TRUE)
}

# Finite-sample item difficulty estimates around a target bank v (≈ \hat b)
# tau is a scale; we can later calibrate it via Fisher info if desired
estimate_bank <- function(v, N = 300, tau = 0.35) {
  v + rnorm(length(v), mean = 0, sd = tau / sqrt(N))
}


# Generate IRT item pool
generate_items <- function(n_items = 100, meanlog = 0.1, sdlog = 0.3, b_sd = 0.75) {
  data.frame(
    a_true = pmax(0.5, rlnorm(n_items, meanlog, sdlog)),
    b_true = rnorm(n_items, mean = 0, sd = b_sd),
    a_cat = NA,
    b_cat = NA
  )
}
# Reusable item pool (same a_true, b_true across sigma/runs)
generate_items_pool <- function(n_items = 100, b_sd = 0.75, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  generate_items(n_items = n_items, b_sd = b_sd)  # your existing generator
}
add_calibration_error <- function(items, sigma,
                                  bias_type  = c("random","structured","none"),
                                  a_bias = 1.1, b_bias = 0.5,
                                  mode       = c("additive","correlated"),
                                  rho = 1.0,
                                  standardize = FALSE
){
  bias_type <- match.arg(bias_type)
  mode      <- match.arg(mode)
  
  b <- items$b_true
  n <- nrow(items)
  
  if (mode == "additive") {
    # OLD BEHAVIOR (backward compatible)
    if (bias_type == "structured") {
      items$b_cat <- b + b_bias + rnorm(length(b), 0, sigma)
      items$a_cat <- items$a_true * a_bias * exp(rnorm(n, 0, sigma/2))
    } else if (bias_type == "random") {
      items$b_cat <- b + rnorm(length(b), 0, sigma)
      items$a_cat <- items$a_true * exp(rnorm(n, 0, sigma))
    } else { # "none"
      items$b_cat <- b
      items$a_cat <- items$a_true
    }
    if (isTRUE(standardize) && bias_type != "none") {
      items$b_cat <- standardize_like(items$b_cat, b)
    }
    
  } else {  # mode == "correlated": control Corr(b, b') = rho, preserve Var(b)
    # Use rho to encode cross-population alignment; ignore sigma by default here.
    mu <- mean(b); s <- sd(b)
    eps <- rnorm(length(b), 0, s)
    b_prime <- mu + rho * (b - mu) + sqrt(pmax(0, 1 - rho^2)) * eps
    if (bias_type == "structured") b_prime <- b_prime + b_bias
    # bias_type == "none": no mean shift; rho handles misalignment
    if (isTRUE(standardize)) b_prime <- standardize_like(b_prime, b)
    items$b_cat <- b_prime
    items$a_cat <- items$a_true
  }
  
  items
}

# Bayesian EAP estimation
estimate_theta <- function(a, b, responses, theta_grid = seq(-4, 4, 0.1)) {
  if (length(responses) < 1) return(0)
  log_likelihood <- sapply(theta_grid, function(theta) {
    prob <- plogis(a * (theta - b))
    sum(responses * log(prob + 1e-10) + (1 - responses) * log(1 - prob + 1e-10))
  })
  likelihood <- exp(log_likelihood - max(log_likelihood))
  prior <- dnorm(theta_grid, 0, 1.5)
  posterior <- likelihood * prior
  if (sum(posterior) < 1e-20) return(0)
  sum(theta_grid * posterior) / sum(posterior)
}


# ---  2. CAT Simulation with Several Stopping Rules ---
# Utility: Compute person-level targeting distance
# Compute targeting distance
compute_targeting_distance <- function(theta_true, item_difficulties) {
  map_dbl(1:length(theta_true), function(i) {
    mean(abs(theta_true[i] - item_difficulties[i, ]), na.rm = TRUE)
  })
}

# --- Full form: administer ALL items to everyone; return final scores like CAT/Fixed ---
run_full_form <- function(
    sigma, b_sd, b_bias = 0, n_persons = 500, n_items = 100,
    bias_type = "random", compute_noisy = FALSE,
    mode = "additive", rho = 0.85, standardize_bprime = FALSE,
    N_cat = 150, N_rescore = 300
){
  items <- generate_items(n_items = n_items, b_sd = b_sd) |>
    add_calibration_error(sigma = sigma, bias_type = bias_type, b_bias = b_bias,
                          mode = mode, rho = rho, standardize = standardize_bprime)
  
  # --- NEW: estimated banks for rescoring ---
  items$bh_true <- estimate_bank(items$b_true, N = N_rescore)
  items$bh_cat  <- estimate_bank(items$b_cat,  N = N_cat)
  
  theta_true <- rnorm(n_persons)
  
  # simulate under TRUE params
  P_true <- plogis(outer(theta_true, items$b_true, "-") *
                     matrix(items$a_true, n_persons, n_items, byrow = TRUE))
  R <- matrix(rbinom(n_persons * n_items, 1, as.vector(P_true)), n_persons, n_items)
  
  theta_hat_true_b  <- numeric(n_persons)
  theta_hat_noisy_b <- if (compute_noisy) numeric(n_persons) else NULL
  
  for (i in 1:n_persons) {
    theta_hat_true_b[i] <- estimate_theta(items$a_true, items$bh_true, R[i, ])    # \hat b
    if (isTRUE(compute_noisy))
      theta_hat_noisy_b[i] <- estimate_theta(items$a_cat, items$bh_cat, R[i, ])   # \hat b'
  }
  
  # targeting vs both banks (for completeness)
  item_difficulties_bprime <- matrix(rep(items$b_cat,  each = n_persons), nrow = n_persons)
  item_difficulties_b      <- matrix(rep(items$b_true, each = n_persons), nrow = n_persons)
  compute_td <- function(mat_b) {
    m <- abs(sweep(mat_b, 1, theta_true, "-"))
    rowMeans(m, na.rm = TRUE)
  }
  
  list(
    theta_true = theta_true,
    theta_hat_true_b  = theta_hat_true_b,     # \hat b
    theta_hat_noisy_b = theta_hat_noisy_b,    # \hat b' (optional)
    test_length = n_items,
    targeting_distance        = compute_td(item_difficulties_b), # primary = TRUE b
    targeting_distance_bprime = compute_td(item_difficulties_bprime),
    targeting_distance_b      = compute_td(item_difficulties_b),
    sigma = sigma, design = "Full", stopping_rule = "full",
    bias_type = bias_type, b_sd = b_sd, b_bias = b_bias
  )
}

# --- NEW: get per-step mean targeting (primary = TRUE bank) -------------------
targeting_means_by_step <- function(theta_true, item_b_by_step) {
  m <- abs(sweep(item_b_by_step, 1, theta_true, "-"))
  colMeans(m, na.rm=TRUE)  # vector: step -> mean distance
}


# --- 2. Main Simulation Wrapper ---
run_combined_simulation <- function(
    # σ knobs
  sigma_draw  = c("range","vector"),
  sigma_range = c(0, 1),
  n_sigma     = 20,
  sigma_vec   = NULL,
  
  # ρ knobs
  rho_draw    = c("vector","range"),
  rho_vec     = NULL,
  rho_range   = c(0, 1),
  n_rho       = 20,
  
  # design axes
  b_sd_vec = c(0.5, 0.75, 1.0),
  stopping_rules = c("SEM","fixed","SEM_topk","full"),
  b_bias_vec = c(-0.5, 0.5),
  n_persons = 100, n_replicates = 30, top_k = 3,
  enforce_max_length = FALSE, full_compute_noisy = FALSE,
  select_bank = "bhat_cat", select_bank_opts = NULL,
  modes = c("correlated"),
  standardize_bprime = FALSE, N_cat = 150, N_rescore = 300,
  
  # NEW: σ handling in correlated
  ignore_sigma_in_correlated = TRUE,
  sigma_anchor_for_correlated = 0,
  
  # I/O
  print_grid = FALSE, return_grid_only = FALSE, grid_seed = NULL
){
  param_grid <- build_param_grid(
    sigma_draw, sigma_range, n_sigma, sigma_vec,
    rho_draw, rho_vec, rho_range, n_rho,
    b_sd_vec, stopping_rules, b_bias_vec, modes,
    select_bank, select_bank_opts,
    ignore_sigma_in_correlated, sigma_anchor_for_correlated,
    grid_seed
  )
  
  if (isTRUE(print_grid)) {
    message("# of base design rows: ", nrow(param_grid))
    print(dplyr::as_tibble(param_grid))
  }
  if (isTRUE(return_grid_only)) return(dplyr::as_tibble(param_grid))
  
  purrr::pmap_dfr(
    param_grid,
    function(sigma, b_sd, stopping_rule, bias_type, b_bias, rho, select_bank, mode) {
      purrr::map_dfr(seq_len(n_replicates), function(seed) {
        set.seed(seed)
        
        res <-
          if (stopping_rule == "fixed") {
            run_fixed_test(
              sigma = sigma, bias_type = bias_type, b_bias = b_bias,
              b_sd = b_sd, n_persons = n_persons,
              mode = mode, rho = ifelse(is.na(rho), 0.85, rho),
              standardize_bprime = standardize_bprime,
              N_cat = N_cat, N_rescore = N_rescore
            )
          } else if (stopping_rule == "full") {
            run_full_form(
              sigma = sigma, b_sd = b_sd, b_bias = b_bias,
              n_persons = n_persons, n_items = 100,
              bias_type = bias_type, compute_noisy = full_compute_noisy,
              mode = mode, rho = ifelse(is.na(rho), 0.85, rho),
              standardize_bprime = standardize_bprime,
              N_cat = N_cat, N_rescore = N_rescore
            )
          } else {
            run_cat_simulation(
              sigma = sigma, stopping_rule = stopping_rule,
              bias_type = bias_type, b_bias = b_bias,
              b_sd = b_sd, n_persons = n_persons,
              top_k = if (stopping_rule == "SEM_topk") top_k else 3,
              enforce_max_length = enforce_max_length,
              select_bank = select_bank,
              mode = mode, rho = ifelse(is.na(rho), 0.85, rho),
              standardize_bprime = standardize_bprime,
              N_cat = N_cat, N_rescore = N_rescore
            )
          }
        
        dplyr::bind_rows(
          calculate_metrics(res, score_with = "noisy", target_on = "auto"),
          calculate_metrics(res, score_with = "true",  target_on = "auto")
        ) |>
          dplyr::mutate(
            enforce_max_length = if (!is.null(res$enforce_max_length)) res$enforce_max_length else NA,
            rho = rho, mode = mode, N_cat = N_cat, N_rescore = N_rescore,
            select_bank = dplyr::coalesce(
              if (!is.null(res$select_bank)) res$select_bank else NA_character_,
              select_bank
            )
          )
      })
    }
  )
}




run_cat_trace <- function(
    sigma, b_sd, b_bias = 0, n_persons = 500, max_steps = 20,
    bias_type = "random", sem_threshold = 0.4, min_steps = 5, top_k = 5
){
  items <- generate_items(100, b_sd = b_sd) |>
    add_calibration_error(sigma = sigma, bias_type = bias_type, b_bias = b_bias)
  
  theta_true <- rnorm(n_persons)
  theta_hat_step_true  <- matrix(NA_real_, n_persons, max_steps)
  theta_hat_step_noisy <- matrix(NA_real_, n_persons, max_steps)
  
  for (i in 1:n_persons) {
    administered <- integer(0); responses <- integer(0); theta_hat <- 0
    for (step in 1:max_steps) {
      avail <- setdiff(seq_len(nrow(items)), administered)
      if (!length(avail)) break
      p  <- plogis(items$a_cat[avail]*(theta_hat - items$b_cat[avail]))
      inf <- items$a_cat[avail]^2 * p * (1 - p)
      pick <- avail[order(inf, decreasing = TRUE)[1:min(top_k, length(inf))]]
      next_item <- if (length(pick)==1) pick else sample(pick, 1, prob = inf[match(pick, avail)])
      
      pr_true <- plogis(items$a_true[next_item]*(theta_true[i] - items$b_true[next_item]))
      responses    <- c(responses, rbinom(1,1,pr_true))
      administered <- c(administered, next_item)
      
      # CAT internal update under b′
      theta_hat <- estimate_theta(items$a_cat[administered], items$b_cat[administered], responses)
      
      # Stepwise scoring under b and under b′
      theta_hat_step_true [i, step] <- estimate_theta(items$a_true[administered], items$b_true[administered], responses)
      theta_hat_step_noisy[i, step] <- estimate_theta(items$a_cat[administered], items$b_cat[administered], responses)
      
      # optional early stop (kept simple; you can add your SEM logic if desired)
      if (step >= min_steps) {
        p_step <- plogis(items$a_cat[administered]*(theta_hat - items$b_cat[administered]))
        info   <- sum(items$a_cat[administered]^2 * p_step * (1 - p_step))
        sem    <- ifelse(info > 0, 1/sqrt(info), Inf)
        if (sem < sem_threshold) break
      }
    }
  }
  
  list(
    theta_true = theta_true,
    theta_hat_step_true  = theta_hat_step_true,
    theta_hat_step_noisy = theta_hat_step_noisy,
    max_steps = max_steps,
    sigma = sigma, b_sd = b_sd, b_bias = b_bias, bias_type = bias_type,
    design = "CAT"
  )
}

# Standard CAT with SEM stopping
#run_cat_simulation(sigma = 1, stopping_rule = "SEM")

# Fixed-length CAT (override early stopping, even SEM rule)
#run_cat_simulation(sigma = 1, stopping_rule = "SEM", enforce_max_length = TRUE)


# --- 3. Fixed-Form Simulation ---
run_fixed_test <- function(
    sigma, b_sd, b_bias = 0, test_length = 20,
    n_persons = 100, bias_type = "random",
    # NEW:
    items_pool = NULL,                   # pass once to keep a_true/b_true fixed across σ
    fixed_indices = NULL,               # pass once to reuse the same fixed form
    mode = c("additive","correlated"),
    rho = 1.0, standardize_bprime = FALSE,
    N_cat = 150, N_rescore = 300
){
  mode <- match.arg(mode)
  items_core <- if (is.null(items_pool)) generate_items_pool(100, b_sd = b_sd) else items_pool
  items <- add_calibration_error(
    items_core,
    sigma        = sigma,
    bias_type    = bias_type,
    b_bias       = b_bias,
    mode         = mode,
    rho          = rho,
    standardize  = standardize_bprime
  )
  items <- estimate_banks(items, N_cat = N_cat, N_rescore = N_rescore)
  
  # Truly fixed form: choose once by a_true and reuse
  if (is.null(fixed_indices)) {
    fixed_indices <- order(-items_core$a_true)[1:test_length]
  }
  
  theta_true <- rnorm(n_persons)
  theta_hat_true_b <- theta_hat_noisy_b <- numeric(n_persons)
  sem_history <- matrix(NA_real_, n_persons, test_length)
  
  item_difficulties_bprime <- matrix(items$b_cat [fixed_indices], n_persons, test_length, byrow = TRUE)
  item_difficulties_b      <- matrix(items$b_true[fixed_indices], n_persons, test_length, byrow = TRUE)
  
  for (i in 1:n_persons) {
    p_true <- plogis(items$a_true[fixed_indices] * (theta_true[i] - items$b_true[fixed_indices]))
    responses <- rbinom(test_length, 1, p_true)
    
    # Rescore with estimated banks
    theta_hat_true_b [i] <- estimate_theta(items$a_true[fixed_indices], items$bh_true[fixed_indices], responses) # \hat b
    theta_hat_noisy_b[i] <- estimate_theta(items$a_cat [fixed_indices], items$bh_cat [fixed_indices], responses) # \hat b′
    
    # optional SEM tracking under the noisy bank
    for (k in 1:test_length) {
      use <- fixed_indices[1:k]
      p_cat <- plogis(items$a_cat[use] * (theta_hat_noisy_b[i] - items$b_cat[use]))
      info  <- sum(items$a_cat[use]^2 * p_cat * (1 - p_cat))
      sem_history[i, k] <- ifelse(info > 0, 1/sqrt(info), Inf)
    }
  }
  
  compute_td <- function(mat_b) {
    m <- abs(sweep(mat_b, 1, theta_true, "-"))
    rowMeans(m, na.rm = TRUE)
  }
  
  list(
    theta_true = theta_true,
    theta_hat_true_b = theta_hat_true_b,
    theta_hat_noisy_b = theta_hat_noisy_b,
    sem_history = sem_history,
    item_difficulties_bprime = item_difficulties_bprime,
    item_difficulties_b = item_difficulties_b,
    targeting_distance        = compute_td(item_difficulties_bprime),
    targeting_distance_bprime = compute_td(item_difficulties_bprime),
    targeting_distance_b      = compute_td(item_difficulties_b),
    test_length = test_length,
    design = "Fixed",
    meta = list(bias_type=bias_type, mode=mode, rho=rho, sigma=sigma, b_bias=b_bias,
                fixed_indices=fixed_indices)
  )
}



# Updated metric calculation
theta_metrics <- function(results, design = "CAT", score_with = c("true", "noisy")) {
  score_with <- match.arg(score_with)
  est <- if (score_with == "true") results$theta_hat_true_b else results$theta_hat_noisy_b
  len <- if (design == "CAT") results$test_lengths else results$test_length
  
  data.frame(
    bias = mean(est - results$theta_true, na.rm = TRUE),
    rmse = sqrt(mean((est - results$theta_true)^2, na.rm = TRUE)),
    cor = cor(est, results$theta_true, use = "complete.obs"),
    avg_length = mean(len),
    sigma = results$sigma,
    b_sd = results$b_sd,
    design = results$design,
    bias_type = results$bias_type,
    stopping_rule = results$stopping_rule,
    enforce_max_length = results$enforce_max_length,
    score_with = score_with,
    targeting_distance = if (!is.null(results$targeting_distance)) mean(results$targeting_distance, na.rm = TRUE) else NA_real_
  )
}


# Standard CAT with SEM stopping
#run_cat_simulation(sigma = 1, stopping_rule = "SEM")

# Fixed-length CAT (override early stopping, even SEM rule)
#run_cat_simulation(sigma = 1, stopping_rule = "SEM", enforce_max_length = TRUE)
run_cat_simulation <- function(
    sigma, b_sd, b_bias = 0, 
    n_persons = 100, max_steps = 20,
    stopping_rule = c("SEM", "fixed", "SEM_topk"),
    sem_threshold = 0.4, min_steps = 5, precision_gain_threshold = 0.02,
    top_k = 3, bias_type = "random", enforce_max_length = FALSE,
    select_bank = c("bhat_cat","bhat_true","bprime","b"),
    # --- NEW ---
    mode = "additive", rho = 0.85, standardize_bprime = FALSE,
    N_cat = 150, N_rescore = 300
){
  stopping_rule <- match.arg(stopping_rule)
  select_bank   <- match.arg(select_bank)
  
  # --- Generate item bank and add calibration error ---
  items <- generate_items(n_items = 100, b_sd = b_sd) |>
    add_calibration_error(sigma = sigma, bias_type = bias_type, b_bias = b_bias,
                          mode = mode, rho = rho, standardize = standardize_bprime)
  
  # --- Estimated banks for rescoring (finite-sample) ---
  items$bh_true <- estimate_bank(items$b_true, N = N_rescore)   # ≈ \hat b
  items$bh_cat  <- estimate_bank(items$b_cat,  N = N_cat)       # ≈ \hat b'
  
  # --- Selection bank (CAT uses this online) ---
  a_sel <- switch(select_bank,
                  "bhat_cat"  = items$a_cat,    # estimated CAT bank
                  "bhat_true" = items$a_true,   # estimated true bank (oracle)
                  "bprime"    = items$a_cat,    # latent misspecified bank (diagnostic)
                  "b"         = items$a_true)   # latent true bank (diagnostic)
  b_sel <- switch(select_bank,
                  "bhat_cat"  = items$bh_cat,
                  "bhat_true" = items$bh_true,
                  "bprime"    = items$b_cat,
                  "b"         = items$b_true)
  
  # --- Simulation setup ---
  theta_true <- rnorm(n_persons)
  theta_hat_true_b  <- numeric(n_persons)   # score with \hat b (rescoring-improved)
  theta_hat_noisy_b <- numeric(n_persons)   # score with \hat b' (admin-coherent)
  test_lengths <- integer(n_persons)
  sem_history  <- matrix(NA_real_, n_persons, max_steps)
  
  # store difficulties per administered step (for diagnostics/targeting)
  item_difficulties_bprime  <- matrix(NA_real_, n_persons, max_steps)
  item_difficulties_b       <- matrix(NA_real_, n_persons, max_steps)
  
  for (i in 1:n_persons) {
    responses <- administered <- NULL
    theta_hat <- 0
    
    for (step in 1:max_steps) {
      available <- setdiff(seq_len(nrow(items)), administered)
      if (!length(available)) break
      p_sel <- plogis(a_sel[available] * (theta_hat - b_sel[available]))
      info  <- a_sel[available]^2 * p_sel * (1 - p_sel)
      next_item <- available[which.max(info)]
      
      # simulate under TRUE params
      p_true <- plogis(items$a_true[next_item] * (theta_true[i] - items$b_true[next_item]))
      responses    <- c(responses, rbinom(1, 1, p_true))
      administered <- c(administered, next_item)
      
      # record TRUE and B' difficulties for targeting diagnostics
      item_difficulties_bprime[i, step] <- items$b_cat [next_item]
      item_difficulties_b     [i, step] <- items$b_true[next_item]
      
      # interim update + SEM under the *selection* bank
      theta_hat <- estimate_theta(a_sel[administered], b_sel[administered], responses)
      p_now     <- plogis(a_sel[administered] * (theta_hat - b_sel[administered]))
      info_now  <- sum(a_sel[administered]^2 * p_now * (1 - p_now))
      current_sem <- ifelse(info_now > 0, 1 / sqrt(info_now), Inf)
      sem_history[i, step] <- current_sem
      
      if (!enforce_max_length && step >= min_steps && stopping_rule == "SEM" && current_sem < sem_threshold) break
      if (enforce_max_length && step == max_steps) break
      if (!enforce_max_length && stopping_rule == "fixed" && step == max_steps) break
    }
    
    test_lengths[i] <- length(administered)
    
    # --- Final post-hoc scoring with estimated banks ---
    theta_hat_true_b [i] <- estimate_theta(items$a_true[administered], items$bh_true[administered], responses)  # rescoring
    theta_hat_noisy_b[i] <- estimate_theta(items$a_cat [administered], items$bh_cat [administered], responses)  # admin-coherent
  }
  
  # --- Targeting diagnostics ---
  compute_td <- function(mat_b) {
    m <- abs(sweep(mat_b, 1, theta_true, "-"))
    rowMeans(m, na.rm = TRUE)
  }
  targeting_distance_bprime <- compute_td(item_difficulties_bprime)
  targeting_distance_b      <- compute_td(item_difficulties_b)
  
  list(
    theta_true = theta_true,
    theta_hat_true_b = theta_hat_true_b,     # ≈ \hat b
    theta_hat_noisy_b = theta_hat_noisy_b,   # ≈ \hat b'
    test_lengths = test_lengths,
    sem_history = sem_history,
    item_difficulties_bprime = item_difficulties_bprime,
    item_difficulties_b = item_difficulties_b,
    targeting_distance        = targeting_distance_bprime,  # legacy
    targeting_distance_bprime = targeting_distance_bprime,
    targeting_distance_b      = targeting_distance_b,
    sigma = sigma, design = "CAT", stopping_rule = stopping_rule,
    bias_type = bias_type, enforce_max_length = enforce_max_length,
    b_sd = b_sd, b_bias = b_bias,
    select_bank = select_bank
  )
}


# --- 3. Fixed-Form Simulation ---
## a. Selects items with highest discrimination parameters
## b. Simulates responses for all examinees at once
## c. Estimates abilities in batch



# --- 4. Metrics + Simulation Runner ---
calculate_metrics <- function(results, design = "CAT",
                              score_with = c("true","noisy"),
                              target_on = "auto") {
  score_with <- match.arg(score_with)
  target_on  <- match.arg(target_on, c("bprime","b","auto"))
  
  est <- if (score_with == "true"  && !is.null(results$theta_hat_true_b))  results$theta_hat_true_b
  else if (score_with == "noisy" && !is.null(results$theta_hat_noisy_b)) results$theta_hat_noisy_b
  else return(data.frame())
  
  len <- if (!is.null(results$test_lengths)) results$test_lengths else results$test_length
  
  td_bprime <- if (!is.null(results$targeting_distance_bprime)) mean(results$targeting_distance_bprime, na.rm=TRUE)
  else if (!is.null(results$targeting_distance))   mean(results$targeting_distance, na.rm=TRUE) else NA_real_
  td_b      <- if (!is.null(results$targeting_distance_b))      mean(results$targeting_distance_b,      na.rm=TRUE) else NA_real_
  td_legacy <- if (target_on == "bprime") td_bprime else if (target_on == "b") td_b else td_b
  
  label <- dplyr::case_when(
    score_with == "noisy" ~ "Admin-coherent (select/score bh_cat)",
    score_with == "true"  ~ "Rescore-improved (select bh_cat, score bh_true)"
  )
  
  dplyr::tibble(
    bias = mean(est - results$theta_true, na.rm=TRUE),
    rmse = sqrt(mean((est - results$theta_true)^2, na.rm=TRUE)),
    pearson_r     = safe_cor(est, results$theta_true, "pearson"),
    spearman_rho  = safe_cor(est, results$theta_true, "spearman"),
    kendall_tau   = safe_cor(est, results$theta_true, "kendall"),
    avg_length = mean(len, na.rm=TRUE),
    sigma = results$sigma, b_sd = results$b_sd,
    design = results$design, bias_type = results$bias_type,
    stopping_rule = results$stopping_rule, score_with = score_with,
    b_bias = if (!is.null(results$b_bias)) results$b_bias else NA_real_,
    select_bank = if (!is.null(results$select_bank)) results$select_bank else NA_character_,
    targeting_distance = td_legacy,
    targeting_distance_bprime = td_bprime,
    targeting_distance_b      = td_b,
    comparison_label = label                                  # <<< NEW
  )
}

################################################ Add two tiny “trace” helpers (stepwise runs) ###############################################
# --- Full-form trace: stepwise scoring under b (and optional b′) ---
run_full_form_trace <- function(
    sigma, b_sd, b_bias = 0,
    n_persons = 400, n_items = 100,
    bias_type = "random",
    # NEW:
    items_pool = NULL, mode = c("additive","correlated"),
    rho = 1.0, standardize_bprime = FALSE
){
  mode <- match.arg(mode)
  items_core <- if (is.null(items_pool)) generate_items_pool(n_items, b_sd = b_sd) else items_pool
  items <- add_calibration_error(
    items_core,
    sigma        = sigma,
    bias_type    = bias_type,
    b_bias       = b_bias,
    mode         = mode,
    rho          = rho,
    standardize  = standardize_bprime
  )
  
  theta_true <- rnorm(n_persons)
  P_true <- plogis(outer(theta_true, items$b_true, "-") *
                     matrix(items$a_true, n_persons, n_items, byrow = TRUE))
  R <- matrix(rbinom(n_persons * n_items, 1, as.vector(P_true)), n_persons, n_items)
  
  theta_hat_step_true <- matrix(NA_real_, n_persons, n_items)
  for (i in 1:n_persons) for (k in 1:n_items) {
    use <- 1:k
    theta_hat_step_true[i, k] <- estimate_theta(items$a_true[use], items$b_true[use], R[i, use])
  }
  list(theta_true = theta_true, theta_hat_step_true = theta_hat_step_true, max_steps = n_items,
       meta = list(bias_type=bias_type, mode=mode, rho=rho, sigma=sigma, b_bias=b_bias))
}
run_fixed_trace <- function(
    sigma, b_sd, b_bias = 0,
    test_length = 20, n_persons = 400,
    bias_type = "random",
    # NEW:
    items_pool = NULL, fixed_indices = NULL,
    mode = c("additive","correlated"),
    rho = 1.0, standardize_bprime = FALSE
){
  mode <- match.arg(mode)
  items_core <- if (is.null(items_pool)) generate_items_pool(100, b_sd = b_sd) else items_pool
  items <- add_calibration_error(
    items_core,
    sigma        = sigma,
    bias_type    = bias_type,
    b_bias       = b_bias,
    mode         = mode,
    rho          = rho,
    standardize  = standardize_bprime
  )
  
  if (is.null(fixed_indices)) {
    fixed_indices <- order(-items_core$a_true)[1:test_length]
  }
  
  theta_true <- rnorm(n_persons)
  theta_hat_step_true <- matrix(NA_real_, n_persons, test_length)
  
  for (i in 1:n_persons) {
    p_true <- plogis(items$a_true[fixed_indices] * (theta_true[i] - items$b_true[fixed_indices]))
    resp   <- rbinom(test_length, 1, p_true)
    
    for (k in 1:test_length) {
      use <- fixed_indices[1:k]
      theta_hat_step_true[i, k] <-
        estimate_theta(items$a_true[use], items$b_true[use], resp[1:k])
    }
  }
  
  list(theta_true = theta_true,
       theta_hat_step_true = theta_hat_step_true,
       max_steps = test_length,
       design = "Fixed",
       meta = list(bias_type=bias_type, mode=mode, rho=rho, sigma=sigma, b_bias=b_bias,
                   fixed_indices=fixed_indices))
}

run_cat_trace_select <- function(
    sigma, b_sd, b_bias = 0, n_persons = 400, max_steps = 30,
    bias_type = "random", sem_threshold = 0.4, min_steps = 5,
    top_k = 5,
    select_bank = c("bhat_cat","bhat_true","bprime","b"),
    score_with = c("noisy","true"),
    # NEW:
    items_pool = NULL,                 # pass a shared pool to keep a_true/b_true fixed across σ
    mode = c("additive","correlated"), # controls how b' is formed
    rho = 1.0,                         # cross-population alignment (1=identical)
    standardize_bprime = FALSE,
    N_cat = 150, N_rescore = 300
){
  select_bank <- match.arg(select_bank)
  score_with  <- match.arg(score_with)
  mode <- match.arg(mode)
  
  # Build or reuse core pool (keeps fixed across runs if provided)
  items_core <- if (is.null(items_pool)) generate_items_pool(100, b_sd = b_sd) else items_pool
  
  # Apply calibration & estimate banks
  items <- add_calibration_error(
    items_core,
    sigma        = sigma,
    bias_type    = bias_type,
    b_bias       = b_bias,
    mode         = mode,
    rho          = rho,
    standardize  = standardize_bprime
  )
  items <- estimate_banks(items, N_cat = N_cat, N_rescore = N_rescore)
  
  # Selection bank
  a_sel <- switch(select_bank,
                  "bhat_cat"  = items$a_cat,
                  "bhat_true" = items$a_true,
                  "bprime"    = items$a_cat,
                  "b"         = items$a_true
  )
  b_sel <- switch(select_bank,
                  "bhat_cat"  = items$bh_cat,   # \hat b′
                  "bhat_true" = items$bh_true,  # \hat b
                  "bprime"    = items$b_cat,    # b′ (oracle)
                  "b"         = items$b_true    # b  (oracle)
  )
  
  # Scoring bank (ops vs. rescoring)
  a_scr <- if (score_with == "noisy") a_sel else items$a_true
  b_scr <- if (score_with == "noisy") b_sel else items$b_true
  
  theta_true           <- rnorm(n_persons)
  theta_hat_step       <- matrix(NA_real_, n_persons, max_steps)
  theta_hat_step_sel   <- matrix(NA_real_, n_persons, max_steps)
  theta_hat_step_true  <- matrix(NA_real_, n_persons, max_steps)
  
  for (i in 1:n_persons) {
    administered <- integer(0); responses <- integer(0); theta_hat <- 0
    for (step in 1:max_steps) {
      avail <- setdiff(seq_len(nrow(items)), administered); if (!length(avail)) break
      p  <- plogis(a_sel[avail] * (theta_hat - b_sel[avail]))
      inf <- a_sel[avail]^2 * p * (1 - p)
      pick <- avail[order(inf, decreasing=TRUE)[1:min(top_k, length(inf))]]
      next_item <- if (length(pick) == 1) pick else sample(pick, 1, prob = inf[match(pick, avail)])
      
      # Response under TRUE DGP
      p_true <- plogis(items$a_true[next_item] * (theta_true[i] - items$b_true[next_item]))
      responses    <- c(responses, rbinom(1, 1, p_true))
      administered <- c(administered, next_item)
      
      # Update selection estimate & record traces
      theta_hat <- estimate_theta(a_sel[administered], b_sel[administered], responses)
      theta_hat_step_sel [i, step] <- theta_hat
      theta_hat_step_true[i, step] <- estimate_theta(items$a_true[administered], items$b_true[administered], responses)
      theta_hat_step     [i, step] <- estimate_theta(a_scr[administered],        b_scr[administered],        responses)
      
      # SEM stop (operational metric)
      if (step >= min_steps) {
        p_now <- plogis(a_sel[administered] * (theta_hat - b_sel[administered]))
        info  <- sum(a_sel[administered]^2 * p_now * (1 - p_now))
        if (info > 0 && 1/sqrt(info) < sem_threshold) break
      }
    }
  }
  
  list(
    theta_true          = theta_true,
    theta_hat_step      = theta_hat_step,       # use this in summarize_rmse_by_step()
    theta_hat_step_sel  = theta_hat_step_sel,
    theta_hat_step_true = theta_hat_step_true,
    max_steps           = max_steps,
    meta = list(bias_type=bias_type, mode=mode, rho=rho, sigma=sigma, b_bias=b_bias)
  )
}


summarize_rmse_by_step <- function(tr){
  steps <- ncol(tr$theta_hat_step_true)
  tibble::tibble(
    step = rep(1:steps, each = length(tr$theta_true)),
    theta_true = rep(tr$theta_true, times = steps),
    est = as.vector(tr$theta_hat_step_true)
  ) |>
    dplyr::group_by(step) |>
    dplyr::summarise(rmse = sqrt(mean((est - theta_true)^2, na.rm = TRUE)), .groups="drop")
}

plot_traj_simple <- function(
    res,
    bias_type    = "none",
    b_bias       = 0,
    n_persons    = 500,
    max_steps    = 25,
    fixed_length = 20,
    methods      = c("cat_ops","cat_rescore","full","fixed_point","fixed_traj"),
    mode         = c("additive","correlated"),
    rho          = NA_real_,
    standardize_bprime = FALSE,
    items_pool   = NULL,
    fixed_indices= NULL,
    N_cat        = 150,
    N_rescore    = 300,
    # --- NEW for stability & clarity ---
    seed         = 20250831,
    n_reps       = 50,
    show_ribbon  = TRUE,
    keep_only    = c("cat_ops","cat_rescore","fixed_traj","full")
){
  stopifnot(is.character(bias_type), length(bias_type)==1)
  mode <- match.arg(mode)
  
  suppressPackageStartupMessages({
    library(dplyr); library(ggplot2); library(tibble); library(purrr)
  })
  
  # ---------- scenario params from `res` ----------
  sub <- res %>% dplyr::filter(.data$bias_type == bias_type, .data$b_bias == b_bias)
  sigma_use <- if (nrow(sub)) stats::median(sub$sigma, na.rm = TRUE) else 0
  b_sd_use  <- if (nrow(sub)) sub$b_sd[which(!is.na(sub$b_sd))[1]] else 0.75
  if (!is.finite(b_sd_use)) b_sd_use <- 0.75
  rho_use   <- if (is.finite(rho)) rho else if ("rho" %in% names(sub) && any(is.finite(sub$rho))) stats::median(sub$rho, na.rm=TRUE) else 0.85
  
  # ---------- shared pool & fixed form ----------
  if (is.null(items_pool))   items_pool   <- generate_items_pool(n_items = 100, b_sd = b_sd_use)
  if (is.null(fixed_indices)) fixed_indices <- order(-items_pool$a_true)[1:fixed_length]
  
  # ---------- labels ----------
  label_map <- c(
    cat_ops     = "CAT (operational)",
    cat_rescore = "CAT (rescoring)",
    full        = "Full-100",
    fixed_point = "Fixed",
    fixed_traj  = sprintf("Fixed-%d", fixed_length)
  )
  
  # ---------- per-rep runner with controlled RNG ----------
  run_one <- function(rep_idx) {
    withr::with_seed(seed + rep_idx, {
      parts <- list()
      
      if ("cat_ops" %in% methods) {
        tr <- run_cat_trace_select(
          sigma = sigma_use, b_sd = b_sd_use, b_bias = b_bias,
          n_persons = n_persons, max_steps = max_steps, bias_type = bias_type,
          select_bank = "bhat_cat", score_with = "noisy",
          items_pool = items_pool, mode = mode, rho = rho_use,
          standardize_bprime = standardize_bprime, N_cat = N_cat, N_rescore = N_rescore
        )
        parts$cat_ops <- summarize_rmse_by_step(tr, rep_id = rep_idx) %>%
          dplyr::mutate(method = label_map[["cat_ops"]])
      }
      
      if ("cat_rescore" %in% methods) {
        tr <- run_cat_trace_select(
          sigma = sigma_use, b_sd = b_sd_use, b_bias = b_bias,
          n_persons = n_persons, max_steps = max_steps, bias_type = bias_type,
          select_bank = "bhat_cat", score_with = "true",
          items_pool = items_pool, mode = mode, rho = rho_use,
          standardize_bprime = standardize_bprime, N_cat = N_cat, N_rescore = N_rescore
        )
        parts$cat_rescore <- summarize_rmse_by_step(tr, rep_id = rep_idx) %>%
          dplyr::mutate(method = label_map[["cat_rescore"]])
      }
      
      if ("full" %in% methods) {
        tr <- run_full_form_trace(
          sigma = sigma_use, b_sd = b_sd_use, b_bias = b_bias,
          n_persons = n_persons, n_items = 100, bias_type = bias_type,
          items_pool = items_pool, mode = mode, rho = rho_use,
          standardize_bprime = standardize_bprime
        )
        parts$full <- summarize_rmse_by_step(tr, rep_id = rep_idx) %>%
          dplyr::mutate(method = label_map[["full"]])
      }
      
      if ("fixed_point" %in% methods) {
        fx <- run_fixed_test(
          sigma = sigma_use, b_sd = b_sd_use, b_bias = b_bias,
          test_length = fixed_length, n_persons = n_persons, bias_type = bias_type,
          items_pool = items_pool, fixed_indices = fixed_indices,
          mode = mode, rho = rho_use, standardize_bprime = standardize_bprime,
          N_cat = N_cat, N_rescore = N_rescore
        )
        rmse_fixed <- sqrt(mean((fx$theta_hat_true_b - fx$theta_true)^2, na.rm = TRUE))
        parts$fixed_point <- tibble::tibble(step = fixed_length, rmse = rmse_fixed,
                                            rep = rep_idx, method = label_map[["fixed_point"]])
      }
      
      if ("fixed_traj" %in% methods) {
        tr <- run_fixed_trace(
          sigma = sigma_use, b_sd = b_sd_use, b_bias = b_bias,
          test_length = fixed_length, n_persons = n_persons, bias_type = bias_type,
          items_pool = items_pool, fixed_indices = fixed_indices,
          mode = mode, rho = rho_use, standardize_bprime = standardize_bprime
        )
        parts$fixed_traj <- summarize_rmse_by_step(tr, rep_id = rep_idx) %>%
          dplyr::mutate(method = label_map[["fixed_traj"]])
      }
      
      dplyr::bind_rows(parts)
    })
  }
  
  # ---------- run reps & keep only desired methods ----------
  df_all <- purrr::map_dfr(seq_len(n_reps), run_one) %>%
    dplyr::filter(.data$method %in% label_map[keep_only])
  
  # ---------- aggregate: mean ± SE over reps ----------
  df_sum <- df_all %>%
    dplyr::group_by(method, step) %>%
    dplyr::summarise(
      mean_rmse = mean(rmse, na.rm = TRUE),
      se        = stats::sd(rmse, na.rm = TRUE) / sqrt(dplyr::n()),
      .groups   = "drop"
    )
  
  # endpoints (last k per method)
  endpts <- df_sum %>%
    dplyr::group_by(method) %>% dplyr::filter(step == max(step)) %>% dplyr::slice_tail(n = 1)
  
  # Full-100 ceiling (dashed line)
  full_row <- df_sum %>% dplyr::filter(method == "Full-100", step == max(step))
  
  title_txt <- if (mode == "correlated") {
    sprintf("Length trajectories under correlated misalignment (ρ = %.2f)", rho_use)
  } else {
    sprintf("Length trajectories under structured shift (Δb = %s)", b_bias)
  }
  
  p <- ggplot(df_sum, aes(step, mean_rmse, color = method, fill = method)) +
    { if (show_ribbon) geom_ribbon(aes(ymin = mean_rmse - se, ymax = mean_rmse + se),
                                   alpha = 0.12, color = NA) } +
    geom_line(linewidth = 1) +
    geom_point(data = endpts, size = 2) +
    scale_x_continuous(breaks = seq(0, max_steps, by = 5), limits = c(0, max_steps)) +
    labs(x = "# items (k)", y = "RMSE",
         title = title_txt,
         subtitle = "CAT-O and CAT-R share selections; divergence reflects scoring only",
         color = NULL, fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  
  if (nrow(full_row)) {
    p <- p +
      geom_hline(yintercept = full_row$mean_rmse[1], linetype = 3) +
      annotate("text", x = max_steps - 1, y = full_row$mean_rmse[1],
               label = "Full-100 ceiling", hjust = 1, vjust = -0.6, size = 3)
  }
  
  p
}

