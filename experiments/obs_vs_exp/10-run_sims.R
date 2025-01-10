library(brlavaan)
library(tidyverse)
library(furrr)
theme_set(theme_bw())
here::i_am("experiments/obs_vs_exp/10-run_sims.R")

## ----- Run simulations -------------------------------------------------------
ncores <- future::availableCores() - 2
future::plan(multisession, workers = ncores)
B <- 250  # Number of simulations

simu_id <-
  expand_grid(
    dist = "Normal",
    model = c("twofac", "growth"),
    rel = 0.8,
    n = c(15, 20, 50, 100, 1000),
    info_penalty = c("observed", "expected"),
    info_bias = c("observed", "expected"),
    info_se = c("observed", "expected")
  ) |>
  rownames_to_column(var = "simid")

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
    lavfun = ifelse(model == "twofac", "sem", "growth"),
    nsimu = B,
    info_penalty = simu_id$info_penalty[i],
    info_bias = simu_id$info_bias[i],
    info_se = simu_id$info_se[i]
  )
  cat("\n")
  save(simu_res, file = "experiments/obs_vs_exp/ove.RData")
}

## ----- Analyse ---------------------------------------------------------------
# load("experiments/obs_vs_exp/ove.RData")
simu_res |>
  imap(\(x, idx) bind_cols(simid = idx, x$simu_res)) |>
  bind_rows()

# # How many converged?
# res |>
#   summarise(
#     fail = any(!converged),
#     .by = c(simu:info_se, method)
#   ) |>
#   drop_na() |>
#   summarise(
#     count = sum(!fail),
#     .by = c(dist:method)
#   ) |>
#   mutate(prop = count / B * 100) |>
#   select(-count) |>
#   pivot_wider(names_from = c(method), values_from = prop) |>
#   print(n = Inf)
#
# # Check bias
# res |>
#   mutate(bias = est - truth) |>
#   summarise(
#     mean_bias = mean(bias, na.rm = TRUE, trim = 0.05),
#     .by = c(dist:info_se, method)
#   ) |>
#   pivot_wider(
#     values_from = mean_bias,
#     names_from = method
#   ) |>
#   print(n = 1000)
#
# # Plot
# res |>
#   mutate(
#     bias = est - truth,
#     across(starts_with("info"), \(x) substr(x, 1, 3)),
#     n = factor(n),
#   ) |>
#   # unite("method", method, starts_with("info"), sep = "-") |>
#   summarise(
#     mean_bias = mean(bias, na.rm = TRUE, trim = 0.1),
#     .by = c(dist, model, rel, n, method, starts_with("info"))
#   ) |>
#   ggplot(aes(as.numeric(n), mean_bias, col = method)) +
#   geom_line() +
#   facet_grid(info_bias * info_se ~ info_penalty)
