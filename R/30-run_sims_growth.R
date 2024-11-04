library(tidyverse)
library(lavaan)
library(furrr)
theme_set(theme_bw())
source("R/10-sem_rbm_functions.R")
source("R/20-gen_data.R")
source("R/21-sim_functions.R")

ncores <- parallel::detectCores() - 2
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
# save(results, file = "R/sim_results.RData")
