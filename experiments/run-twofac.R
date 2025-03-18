source(here::here("experiments/_setup.R"))

simu_res_twofac <- vector("list", length = nrow(simu_id))
for (i in seq_len(nrow(simu_id))) {
  dist  <- simu_id$dist[i]
  model <- "twofac"
  rel   <- simu_id$rel[i]
  n     <- simu_id$n[i]
  seeds <- simu_id$seed[[i]][1:B]
  # seeds <- NULL

  cli::cli_inform(">>> {Sys.time()} <<<\n\n[{i} / {nrow(simu_id)}] Now running {model} models ({dist}) rel = {rel}, n = {n}\n")
  simu_res_twofac[[i]] <- sim_fun(
    dist = dist,
    model = model,
    rel = rel,
    n = n,
    nsimu = B,
    lavsim = FALSE,
    whichsims = c("ML", "eRBM", "iRBM"),
    bounds = "standard",
    keep_going = FALSE,
    data_scale = 1,
    seeds = seeds,
    maxgrad = FALSE
  )
  cat("\n")
  save(simu_res_twofac, file = "experiments/simu_res_twofac.RData")
}

# map(simu_res_twofac, \(x) x$simu_res) |>
#   bind_rows() |>
#   summarise(
#     count = sum(converged) / B * 100,
#     .by = c(dist:method)
#   ) |>
#   pivot_wider(names_from = method, values_from = count) |>
#   print(n = Inf)
