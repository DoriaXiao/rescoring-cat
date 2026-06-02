# ================== Side-by-side: Additive vs Correlated ==================
library(dplyr); library(ggplot2); library(purrr); library(tidyr)
source("R_ml/simulation_ml_source.R")
tr1 <- run_cat_trace_select(sigma = 0.4, b_sd = 0.75, bias_type = "random",
                            select_bank = "bhat_cat", score_with = "true",
                            n_persons = 400, max_steps = 25, do_stop = FALSE)
tr2 <- run_cat_trace_select(sigma = 0.4, b_sd = 0.75, bias_type = "structured",
                            select_bank = "bhat_true", score_with = "true",
                            n_persons = 400, max_steps = 25, do_stop = FALSE)

train_df <- dplyr::bind_rows(
  make_mlstop_training(tr1, TRUE, 0.01),
  make_mlstop_training(tr2, TRUE, 0.01)
)

train_df$bias_type   <- factor(train_df$bias_type)
train_df$select_bank <- factor(train_df$select_bank)

ml_stop_model <- glm(
  y ~ sem + dsem + step + td_bprime + td_b + sigma + rho + b_sd +
    bias_type + select_bank,
  data = train_df,
  family = binomial()
)
summary(ml_stop_model)

# Optionally persist for reuse
saveRDS(ml_stop_model, file = "ml_stop_model.rds")


ml_stop_model <- readRDS("ml_stop_model.rds")

res <- run_cat_simulation(
  sigma = 0.4, b_sd = 0.75, bias_type = "random",
  stopping_rule = "ML_stop",
  ml_stop_model = ml_stop_model, ml_p_cut = 0.5,
  select_bank = "bhat_cat",
  n_persons = 500, max_steps = 25, rho = 0.85
)

# Evaluate
dplyr::bind_rows(
  calculate_metrics(res, score_with = "noisy"),
  calculate_metrics(res, score_with = "true")
)


compare_rules <- function(...) {
  rules <- c("SEM", "fixed", "ML_stop")
  dplyr::bind_rows(lapply(rules, function(r){
    res <- run_cat_simulation(
      stopping_rule = r,
      ml_stop_model = if (r=="ML_stop") ml_stop_model else NULL,
      ...
    )
    dplyr::bind_rows(
      calculate_metrics(res, score_with = "true"),
      calculate_metrics(res, score_with = "noisy")
    ) |> dplyr::mutate(rule = r)
  }))
}

cmp <- compare_rules(
  sigma = 0.4, b_sd = 0.75, bias_type = "random",
  select_bank = "bhat_cat", n_persons = 500
)

library(ggplot2)
ggplot(cmp, aes(avg_length, rmse, color = rule, shape = score_with)) +
  geom_point(size = 3) +
  labs(x = "Average Test Length", y = "RMSE (θ̂ vs θ)")


##############################
# a) Enrich training data

sigs <- c(0.2, 0.4, 0.8)
rhos <- c(0.6, 0.85)
b_sds <- c(0.5, 0.75, 1.0)


# b) stronger learner
library(randomForest)
train_df2 <- train_df
train_df2$dsem[is.na(train_df2$dsem)] <- 0
train_df2$td_bprime[is.na(train_df2$td_bprime)] <- mean(train_df2$td_bprime, na.rm = TRUE)
train_df2$td_b[is.na(train_df2$td_b)]           <- mean(train_df2$td_b, na.rm = TRUE)

rf_model <- randomForest(
  factor(y) ~ sem + dsem + step + td_bprime + td_b + sigma + rho + b_sd,
  data = train_df2, ntree = 300
)

print(rf_model)
varImpPlot(rf_model)

saveRDS(rf_model, "rf_stop_model.rds")
rf_model <- readRDS("rf_stop_model.rds")  # later when running sims
cmp <- compare_rules(
  ml_model = rf_model,          # <- pass here
  ml_p_cut = 0.3,
  sigma = 0.4, b_sd = 0.75,     # the rest go in ...
  bias_type = "random",
  select_bank = "bhat_cat",
  n_persons = 500
)

ggplot(cmp, aes(avg_length, rmse, color = rule, shape = score_with)) +
  geom_point(size = 3) +
  labs(x = "Average Test Length", y = "RMSE (θ̂ vs θ)")


# Example use (same design you’ve been plotting)
ml_curve <- sweep_ml_cuts(
  cuts = seq(0.2, 0.9, by = 0.05),
  rf_model = rf_model,
  sigma = 0.4, b_sd = 0.75, bias_type = "random",
  select_bank = "bhat_cat", n_persons = 500, max_steps = 25, rho = 0.85
)

library(ggplot2)
ggplot(ml_curve, aes(ml_p_cut, avg_length, shape = score_with)) +
  geom_point() + geom_line() +
  labs(x = "ml_p_cut", y = "Average Test Length (ML_stop)")


# First, measure SEM at this design
sem_res <- run_cat_simulation(
  stopping_rule = "SEM",
  sigma = 0.4, b_sd = 0.75, bias_type = "random",
  select_bank = "bhat_cat", n_persons = 500, max_steps = 25, rho = 0.85
)
sem_len <- calculate_metrics(sem_res, score_with = "true")$avg_length

# Find ml_p_cut whose avg length is closest to SEM’s
find_cut_for_length <- function(target_len, cuts, rf_model, ...) {
  grid <- sweep_ml_cuts(cuts, rf_model, ...)
  grid_true <- subset(grid, score_with == "true")
  grid_true$len_gap <- abs(grid_true$avg_length - target_len)
  grid_true[which.min(grid_true$len_gap), c("ml_p_cut","avg_length","rmse")]
}

best <- find_cut_for_length(
  target_len = sem_len,
  cuts = seq(0.2, 0.9, by = 0.05),
  rf_model = rf_model,
  sigma = 0.4, b_sd = 0.75, bias_type = "random",
  select_bank = "bhat_cat", n_persons = 500, max_steps = 25, rho = 0.85
)
best


# SEM baseline
sem_res <- run_cat_simulation(
  stopping_rule = "SEM",
  sigma = 0.4, b_sd = 0.75, bias_type = "random",
  select_bank = "bhat_cat", n_persons = 500, max_steps = 25, rho = 0.85
)

# ML_stop with your chosen cutoff
ml_res <- run_cat_simulation(
  stopping_rule = "ML_stop",
  ml_stop_model = rf_model,
  ml_p_cut = 0.30,
  sigma = 0.4, b_sd = 0.75, bias_type = "random",
  select_bank = "bhat_cat", n_persons = 500, max_steps = 25, rho = 0.85
)

dplyr::bind_rows(
  calculate_metrics(sem_res, score_with = "true") |> dplyr::mutate(rule = "SEM"),
  calculate_metrics(ml_res,  score_with = "true") |> dplyr::mutate(rule = "ML_stop (0.30)")
)[, c("rule","avg_length","rmse","bias","pearson_r")]
