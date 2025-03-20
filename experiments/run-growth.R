source(here::here("experiments/_setup.R"))

simu_res_growth <- vector("list", length = nrow(simu_id))
for (i in seq_len(nrow(simu_id))) {
  dist  <- simu_id$dist[i]
  model <- "growth"
  rel   <- simu_id$rel[i]
  n     <- simu_id$n[i]
  seeds <- simu_id$seed[[i]][1:B]
  # seeds <- NULL

  cli::cli_inform(">>> {Sys.time()} <<<\n\n[{i} / {nrow(simu_id)}] Now running {model} models ({dist}) rel = {rel}, n = {n}\n")
  simu_res_growth[[i]] <- sim_fun(
    dist = dist,
    model = model,
    rel = rel,
    n = n,
    nsimu = B,
    lavsim = FALSE,
    whichsims = c("ML", "eRBM", "iRBM"),
    bounds = "standard",
    keep_going = FALSE,
    data_scale = 1 / 10,
    seeds = seeds,
    maxgrad = FALSE
  )
  cat("\n")
  save(simu_res_growth, file = "experiments/simu_res_growth_new.RData")
}

# map(simu_res_growth, \(x) x$simu_res) |>
#   bind_rows() |>
#   summarise(
#     count = sum(converged, na.rm = TRUE) / B * 100,
#     .by = c(dist:method)
#   ) |>
#   pivot_wider(names_from = method, values_from = count) |>
#   print(n = Inf)
#
# simu_res_growth[[1]]$simu_res |>
#   mutate(param = map(est, names), .before = est) |>
#   unnest(c(param, est, se, truth)) |>
#   drop_na(param) |>
#   filter(converged, !is.na(se)) |>
#   summarise(
#     bias = mean((est - truth) / truth, na.rm = TRUE),
#     .by = dist:param
#   ) |>
#   mutate(method = factor(method, levels = c("ML", "eRBM", "iRBM"))) |>
#   ggplot(aes(param, bias, fill = method)) +
#   geom_col(position = "dodge")
