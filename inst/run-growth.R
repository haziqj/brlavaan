source(here::here("experiments/_setup.R"))

simu_res_growth <- vector("list", length = nrow(simu_id))
for (i in c(1, 6)) {
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
    bounds = "standard",
    data_scale = 1 / 10,
    seeds = seeds
  )
  cat("\n")
  # save(simu_res_growth, file = "experiments/simu_res_growth.RData")
}

bind_rows(simu_res_growth) |>
  mutate(seOK = map_lgl(se, ~!any(is.na(.x)))) |>
  summarise(
    count = sum(converged & seOK & Sigma_OK, na.rm = TRUE) / B * 100,
    .by = c(dist:method)
  ) |>
  pivot_wider(names_from = method, values_from = count) |>
  print(n = Inf)

mycols <- c(
  ML = "#E31A1C",
  lav = "#FB9A99",
  eRBM = "#A6CEE3",
  iRBM = "#1F78B4",
  JB = "#B2DF8A",
  BB = "#33A02C",
  REML = "#FDBF6F",
  Ozenne = "#FF7F00"
)

bind_rows(
  list(
    "none" = bind_rows(simu_res_growth),
    "standard" = stdbounds
  ),
  .id = "bounds"
) |>
  mutate(param = map(truth, names), .before = est) |>
  unnest(c(param, est, se, truth)) |>
  filter(converged, !is.na(se)) |>
  summarise(
    bias = mean((est - truth) / truth, na.rm = TRUE, trim = 0),
    .by = c(dist:param, bounds)
  ) |>
  mutate(
    method = factor(method, levels = names(mycols)),
    param = factor(param, levels = c("i~~i", "i~1", "s~~s", "s~1", "i~~s", "v")),
    bounds = factor(bounds, levels = c("none", "standard"), labels = c("bounds = \"none\"", "bounds = \"standard\"")),
    rel = factor(rel, levels = c("0.8", "0.5"), labels = paste0("Rel = ", c("0.8", "0.5"))),
  ) |>
  filter(param %in% c("i~~i", "s~~s",  "i~~s", "v")) |>
  ggplot(aes(param, abs(bias), fill = method)) +
  geom_col(position = "dodge", width = 0.7) +
  # geom_hline(yintercept = 0, linetype = 2) +
  facet_grid(rel ~ bounds) +
  labs(x = NULL, y = "Absolute Relative Bias", fill = NULL, title = "GROWTH CURVE MODEL", subtitle = "Normal, n = 15, repl. = 500", caption = "ML uses fit_sem(); lav uses growth(); JB = Jacknife; BB = Bootstrap.") +
  scale_fill_manual(values = mycols) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_x_discrete(labels = c(
    expression(Psi["1,1"]),
    # expression(alpha[1]),
    expression(Psi["2,2"]),
    # expression(alpha[2]),
    expression(Psi["1,2"]),
    expression(Theta["1,1"])
  ))
