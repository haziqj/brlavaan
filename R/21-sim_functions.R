library(tidyverse)
library(furrr)
library(lavaan)
source("R/10-sem_rbm_functions.R")
source("R/20-gen_data.R")

ncores <- parallel::detectCores() - 2
future::plan(multisession, workers = ncores)
nsimu <- 4

sim_fun <- function(dist = "Normal", model = "twofac", rel = 0.8, n = 10, lavsim = FALSE) {
  dist <- match.arg(dist, c("Normal", "Kurtosis", "Non-normal"))
  model <- match.arg(model, c("growth", "twofac"))
  rel <- match.arg(as.character(rel), c("0.8", "0.5"))

  # Generate all data ----------------------------------------------------------
  if (model == "growth") {
    gen_data <- gen_data_growth
    txt_mod <- txt_mod_growth
  } else {
    gen_data <- gen_data_twofac
    txt_mod <- txt_mod_twofac
  }
  datasets <- replicate(
    nsimu,
    gen_data(n = n, rel = rel, dist = dist, lavsim = lavsim),
    simplify = FALSE
  )
  mod <- txt_mod(rel)
  true_vals <- truth(datasets[[1]])

  simu_res <- future_map(
    seq_len(nsimu),
    safely(function(j) {
      dat <- datasets[[j]]

      fit_ML    <- fit_sem(model = mod, data = dat, method = "ML")
      fit_eRBM  <- fit_sem(model = mod, data = dat, method = "eRBM")
      fit_iRBM  <- fit_sem(model = mod, data = dat, method = "iRBM")
      fit_iRBMp <- fit_sem(model = mod, data = dat, method = "iRBMp")

      data.frame(
        dist = dist,
        model = model,
        rel = rel,
        n = n,
        param = rep(names(true_vals), 4),
        est = c(fit_ML$coef, fit_eRBM$coef, fit_iRBM$coef, fit_iRBMp$coef),
        se = c(fit_ML$stderr, fit_eRBM$stderr, fit_iRBM$stderr, fit_iRBMp$stderr),
        truth = rep(true_vals, 4),
        method = rep(c("ML", "eRBM", "iRBM", "iRBMp"), each = length(true_vals)),
        time = rep(c(fit_ML$time, fit_eRBM$time, fit_iRBM$time, fit_iRBMp$time), each = length(true_vals))
      )
    }, otherwise = NA),
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
  )
  simu_res

  # Clean up -------------------------------------------------------------------
  where_error <- which(sapply(simu_res, \(x) !is.null(x$error)))
  error_list <- lapply(simu_res[where_error], \(x) x$error)
  names(error_list) <- where_error

  list(
    simu_res = do.call("rbind", lapply(simu_res, \(x) x$result)),
    errors = error_list
  )
}

# simu_id <-
#   expand_grid(
#     dist = c("Normal", "Kurtosis", "Non-normal"),
#     model = c("growth", "twofac"),
#     rel = factor(c(0.8, 0.5), levels = c(0.8, 0.5)),
#     n = c(15, 20, 50, 100, 1000),
#   ) |>
#   rownames_to_column(var = "simid")

simu_id <-
  expand_grid(
    # dist = c("Normal", "Kurtosis", "Non-normal"),
    dist = "Normal",
    model = "twofac",
    rel = c(0.8, 0.5),
    n = c(15, 20, 50, 100, 1000),
  ) |>
  rownames_to_column(var = "simid")

simu_res <- vector("list", length = nrow(simu_id))
for (i in seq_len(nrow(simu_id))) {
  dist  <- simu_id$dist[i]
  model <- simu_id$model[i]
  rel   <- simu_id$rel[i]
  n     <- simu_id$n[i]

  cli::cli_inform("[{i} / {nrow(simu_id)}] Now running {model} models ({dist}) rel = {rel}, n = {n}\n")
  simu_res[[i]] <- sim_fun(dist, model, rel, n, lavsim = TRUE)
  cat("\n")
}

results <- do.call(rbind, lapply(simu_res, \(x) x$simu_res))
results |>
  group_by(dist, model, rel, n, method) |>
  # filter(abs(est) < 300) |>
  summarise(
    bias = mean(est - truth),
    mse = mean((est - truth) ^ 2),
    cover = mean((truth - 1.96 * se < est) & (est < truth + 1.96 * se), na.rm = TRUE),
    time = mean(time),
    nsimu = n()
  ) |>
  filter(rel == 0.8) |>
  ggplot(aes(x = as.numeric(factor(n)), y = (bias), col = method)) +
  geom_line() +
  facet_grid( ~ dist)

