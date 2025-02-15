simu_res_twofac <- vector("list", length = nrow(simu_id))
for (i in seq_len(nrow(simu_id))) {
  dist  <- simu_id$dist[i]
  model <- simu_id$model[i]
  rel   <- simu_id$rel[i]
  n     <- simu_id$n[i]
  seeds <- simu_id$seed[i]

  cli::cli_inform(">>> {Sys.time()} <<<\n\n[{i} / {nrow(simu_id)}] Now running {model} models ({dist}) rel = {rel}, n = {n}\n")
  simu_res_twofac[[i]] <- sim_fun(
    dist = dist,
    model = model,
    rel = rel,
    n = n,
    nsimu = B,
    lavsim = FALSE,
    whichsims = c("ML", "eRBM", "iRBM"),
    info_pen = "observed",
    info_bias = "observed",
    info_se = "observed",
    seeds = seeds,
    keep_going = FALSE,
    data_scale = 1
  )
  cat("\n")
  save(simu_res_twofac, file = "experiments/simu_res_twofac.RData")
}

