# 1/7/25: Want to rerun experiments with sandwich vcov for improved SEs

source(here::here("experiments/brsem/_setup.R"))

simu_res_serobust_growth <- vector("list", length = nrow(simu_id))
for (i in seq_len(nrow(simu_id))) {
  dist  <- simu_id$dist[i]
  model <- "growth"
  rel   <- simu_id$rel[i]
  n     <- simu_id$n[i]
  seeds <- simu_id$seed[[i]][1:B]
  # seeds <- NULL

  cli::cli_inform(">>> {Sys.time()} <<<\n\n[{i} / {nrow(simu_id)}] Now running {model} models ({dist}) rel = {rel}, n = {n}\n")
  simu_res_serobust_growth[[i]] <- sim_fun(
    dist = dist,
    model = model,
    rel = rel,
    n = n,
    nsimu = B,
    nboot = 500L,
    lavsim = FALSE,
    whichsims = c("ML", "eRBM", "iRBM", "Ozenne"),
    bounds = "none",
    information = "observed",
    se = "robust.huber.white",
    data_scale = 1 / 10,
    seeds = seeds
  )
  cat("\n")
  save(simu_res_serobust_growth, file = here::here("experiments/brsem/simu_res_serobust_growth.RData"))
}

# bind_rows(simu_res_serobust_growth) |>
#   mutate(seOK = map_lgl(se, ~!any(is.na(.x)))) |>
#   summarise(
#     count = sum(converged & seOK & Sigma_OK, na.rm = TRUE) / B * 100,
#     .by = c(dist:method)
#   ) |>
#   pivot_wider(names_from = method, values_from = count) |>
#   print(n = Inf)
#
#
# simu_res_serobust_growth[[1]] |>
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
