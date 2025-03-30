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
    nboot = 500L,
    lavsim = FALSE,
    whichsims = c("ML", "eRBM", "iRBM", "Ozenne", "REML", "JB", "BB"),
    bounds = "none",
    data_scale = 1,
    seeds = seeds
  )
  cat("\n")
  save(simu_res_twofac, file = here::here("experiments/simu_res_twofac.RData"))
}

# bind_rows(simu_res_growth) |>
#   mutate(seOK = map_lgl(se, ~!any(is.na(.x)))) |>
#   summarise(
#     count = sum(converged & seOK & Sigma_OK, na.rm = TRUE) / B * 100,
#     .by = c(dist:method)
#   ) |>
#   pivot_wider(names_from = method, values_from = count) |>
#   print(n = Inf)
#
#
# simu_res_twofac[[1]] |>
#   mutate(param = map(est, names), .before = est) |>
#   unnest(c(param, est, se, truth)) |>
#   drop_na(param) |>
#   filter(abs(est) < 10) |>
#   summarise(
#     bias = mean((est - truth) / truth, na.rm = TRUE),
#     .by = dist:param
#   ) |>
#   mutate(method = factor(method, levels = c("ML", "eRBM", "iRBM"))) |>
#   ggplot(aes(param, bias, fill = method)) +
#   geom_col(position = "dodge")
