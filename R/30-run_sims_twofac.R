library(tidyverse)
library(lavaan)
library(furrr)
theme_set(theme_bw())
source("R/12-sem_rbm_functions.R")
source("R/20-gen_data.R")
source("R/21-sim_functions.R")

ncores <- parallel::detectCores() - 2
future::plan(multisession, workers = ncores)
nsimu <- 150

simu_id <-
  expand_grid(
    dist = c("Normal", "Kurtosis", "Non-normal"),
    model = c("twofac"),
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

  cli::cli_inform(">>> {Sys.time()} <<<\n\n[{i} / {nrow(simu_id)}] Now running {model} models ({dist}) rel = {rel}, n = {n}\n")
  simu_res[[i]] <- sim_fun(
    dist = dist,
    model = model,
    rel = rel,
    n = n,
    lavsim = FALSE,
    lavfun = "cfa",
    meanstructure = FALSE,
    bounds = "standard"
  )
  cat("\n")
}

res_twofac <-
  do.call(rbind, lapply(simu_res, \(x) x$simu_res)) |>
  mutate(
    method = factor(method, levels = c("ML", "eRBM", "iRBM", "iRBMp")),
    dist = factor(dist, levels = c("Normal", "Kurtosis", "Non-normal")),
    # model = factor(model, labels = c("Growth model", "Two factor model")),
    rel = factor(rel, levels = c(0.8, 0.5), labels = c("Rel = 0.8", "Rel = 0.5")),
    n = factor(n)
  )
save(res_twofac, file = "R/sim_results_twofac.RData")

res_twofac |>
  as_tibble() |>
  dplyr::filter(dist == "Normal") |>
  dplyr::filter(param %in% c("y1~~y1", "fx~~fx", "fy~~fy", "fy~fx", "fx=~x2")) |>
  summarise(
    truth = first(truth),
    bias = mean(est - truth, na.rm = TRUE),
    .by = c(rel, n, param, method)
  ) |>
  mutate(bias = bias / truth) |>
  pivot_wider(names_from = method, values_from = bias)
