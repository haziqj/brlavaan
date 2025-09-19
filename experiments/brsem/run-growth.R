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
    nboot = 500L,
    lavsim = FALSE,
    whichsims = c("ML", "eRBM", "iRBM", "Ozenne", "REML", "JB", "BB"),
    bounds = "none",
    data_scale = 1 / 10,
    seeds = seeds
  )
  cat("\n")
  save(simu_res_growth, file = here::here("experiments/brsem/simu_res_growth.RData"))
}
#
# bind_rows(simu_res_growth) |>
#   mutate(seOK = map_lgl(se, ~!any(is.na(.x)))) |>
#   summarise(
#     count = sum(converged & seOK & Sigma_OK, na.rm = TRUE) / B * 100,
#     .by = c(dist:method)
#   ) |>
#   pivot_wider(names_from = method, values_from = count) |>
#   print(n = Inf)
#
# mycols <- RColorBrewer::brewer.pal(8, "Paired")
# names(mycols) <- c("eRBM", "iRBM", "JB", "BB", "lav", "ML", "REML", "Ozenne")
#
# bind_rows(simu_res_growth) |>
#   unnest(c(est, se, truth)) |>
#   mutate(param = names(truth)) |>
#   distinct(seed, sim, dist, model, rel, n, method, param, .keep_all = TRUE) |>
#   filter(param %in% c("i~~i", "s~~s", "i~~s", "v"), !is.na(se), Sigma_OK, converged) |>
#   summarise(
#     bias = mean((est - truth) / truth, trim = 0.05),
#     .by = c(dist, model, rel, n, method, param)
#   ) |>
#   ggplot(aes(n, bias, col = method)) +
#   geom_line() +
#   facet_grid(param ~ rel + dist) +
#   scale_x_log10() +
#   scale_colour_manual(values = mycols)
#
