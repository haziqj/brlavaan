simu_id <-
  expand_grid(
    dist = c("Normal", "Kurtosis", "Non-normal"),
    model = "growth",
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

  cli::cli_inform(">>> {Sys.time()} <<<\n\n[{i} / {nrow(simu_id)}] Now running {model} models ({dist}) rel = {rel}, n = {n}\n")
  simu_res[[i]] <- sim_fun(
    dist = dist,
    model = model,
    rel = rel,
    n = n,
    lavsim = FALSE,
    lavfun = "growth",
    nsimu = B
  )
  cat("\n")
}

# res_growth <-
#   do.call(rbind, lapply(simu_res, \(x) x$simu_res)) |>
#   mutate(
#     method = factor(method, levels = c("ML", "eRBM", "iRBM", "iRBMp")),
#     dist = factor(dist, levels = c("Normal", "Kurtosis", "Non-normal")),
#     model = factor(model, labels = c("Growth model", "Two factor model")),
#     rel = factor(rel, levels = c(0.8, 0.5), labels = c("Rel = 0.8", "Rel = 0.5")),
#     n = factor(n)
#   )
simu_res_growth <- simu_res
save(simu_res_growth, file = "experiments/simu_res_growth.RData")
