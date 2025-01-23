simu_id <-
  expand_grid(
    dist = c("Normal", "Kurtosis", "Non-normal"),
    model = "twofac",
    rel = c(0.8, 0.5),
    n = c(15, 20, 50, 100, 1000)
    # dist = "Non-normal",
    # model = "twofac",
    # rel = 0.5,
    # n = 20
  ) |>
  rownames_to_column(var = "simid")

simu_res_twofac <- vector("list", length = nrow(simu_id))
for (i in seq_len(nrow(simu_id))) {
  dist  <- simu_id$dist[i]
  model <- simu_id$model[i]
  rel   <- simu_id$rel[i]
  n     <- simu_id$n[i]

  cli::cli_inform(">>> {Sys.time()} <<<\n\n[{i} / {nrow(simu_id)}] Now running {model} models ({dist}) rel = {rel}, n = {n}\n")
  simu_res_twofac[[i]] <- sim_fun(
    dist = dist,
    model = model,
    rel = rel,
    n = n,
    nsimu = B,
    lavsim = FALSE,
    lavfun = "growth",
    whichsims = c("ML", "eRBM", "iRBM"),
    info_pen = "observed",
    info_bias = "observed",
    info_se = "observed"
  )
  cat("\n")
  save(simu_res_twofac, file = "experiments/simu_res_twofac.RData")
}

