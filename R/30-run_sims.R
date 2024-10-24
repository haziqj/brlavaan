library(tidyverse)
library(lavaan)
library(furrr)
theme_set(theme_bw())
source("R/10-sem_rbm_functions.R")
source("R/20-gen_data.R")
source("R/21-sim_functions.R")

ncores <- parallel::detectCores() - 1
future::plan(multisession, workers = ncores)
nsimu <- 1000

simu_id <-
  expand_grid(
    dist = c("Normal", "Kurtosis", "Non-normal"),
    model = c("twofac", "growth"),
    rel = c(0.8, 0.5),
    n = c(15, 20, 50, 100, 1000),
  ) |>
  rownames_to_column(var = "simid")

# Run simulations --------------------------------------------------------------
simu_res <- vector("list", length = nrow(simu_id))
for (i in seq_len(nrow(simu_id))) {
  dist  <- simu_id$dist[i]
  model <- simu_id$model[i]
  rel   <- simu_id$rel[i]
  n     <- simu_id$n[i]

  cli::cli_inform("[{i} / {nrow(simu_id)}] Now running {model} models ({dist}) rel = {rel}, n = {n}\n")
  simu_res[[i]] <- sim_fun(dist, model, rel, n, lavsim = FALSE)
  cat("\n")
}

results <-
  do.call(rbind, lapply(simu_res, \(x) x$simu_res)) |>
  mutate(
    method = factor(method, levels = c("ML", "eRBM", "iRBM", "iRBMp")),
    dist = factor(dist, levels = c("Normal", "Kurtosis", "Non-normal")),
    model = factor(model, labels = c("Growth model", "Two factor model")),
    rel = factor(rel, levels = c(0.8, 0.5), labels = c("Rel = 0.8", "Rel = 0.5")),
    n = factor(n)
  )
save(results, file = "R/sim_results.RData")

# Analysis ---------------------------------------------------------------------
load("R/sim_results.RData")

# Overall plot
results |>
  drop_na() |>
  filter(abs(est - truth) / truth < 2) |>
  # filter(method %in% c("iRBM", "ML")) |>
  group_by(dist, model, rel, n, method) |>
  summarise(bias = mean(est - truth, na.rm = TRUE, trim = 0.01)) |>
  ggplot(aes(x = as.numeric(n), y = abs(bias), col = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "gray30") +
  geom_line(size = 0.8) +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000)) +
  facet_grid(model + rel ~ dist, scales = "free_y") +
  labs(
    x = "Sample size (n)",
    y = "Absolute mean bias",
    col = NULL
  ) +
  theme(legend.position = "top")

# Two factor model results
results |>
  filter(model == "Two factor model", dist %in% c("Normal", "Non-normal")) |>
  filter(param %in% c("y1~~y1", "fx~~fx", "fy~~fy", "fy~fx", "fx=~x2")) |>
  filter(abs(est) < 5) |>
  filter(method != "iRBMp") |>
  filter(method != "iRBM") |>

  group_by(param, dist, rel, n, method) |>
  summarise(
    bias = mean(est - truth, na.rm = TRUE, trim = 0.01),
    truth = first(truth)
  ) |>
  mutate(
    rel_bias = bias / truth,
    param = factor(param, levels = c("y1~~y1", "fx~~fx", "fy~~fy", "fy~fx", "fx=~x2")),
  ) |>
  ggplot(aes(x = as.numeric(n), y = rel_bias, col = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "gray30") +
  geom_line() +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000)) +
  facet_grid(param ~ rel + dist, scales = "free_y") +
  labs(
    x = "Sample size (n)",
    y = "Relative mean bias",
    col = NULL
  ) +
  theme(legend.position = "top")

results |>
  filter(model == "Two factor model") |>
  filter(param %in% c("y1~~y1", "fx~~fx", "fy~~fy", "fy~fx", "fx=~x2")) |>
  filter(abs(est) < 5) |>
  group_by(dist, rel, param, n, method) |>
  summarise(
    # bias = mean(est - truth, na.rm = TRUE),
    mse = mean((est - truth) ^ 2, na.rm = TRUE),
    # cover = mean((truth - 1.96 * se < est) & (est < truth + 1.96 * se), na.rm = TRUE),
    # time = mean(time, na.rm = TRUE),
    # nsimu = n()
  ) |>
  pivot_wider(names_from = n, values_from = mse, names_prefix = "n = ") |>
  print(n = 100)

