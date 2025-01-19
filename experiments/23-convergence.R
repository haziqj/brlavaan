simu_res <- c(simu_res_growth, simu_res_twofac)
simu_id <-
  expand_grid(
    model = c("growth", "twofac"),
    dist = c("Normal", "Kurtosis", "Non-normal"),
    rel = c(0.8, 0.5),
    n = c(15, 20, 50, 100, 1000)
  ) |>
  rownames_to_column(var = "simid")

# Create results data frame
res <-
  simu_res |>
  imap(\(x, idx) bind_cols(simid = idx, x$simu_res)) |>
  bind_rows() |>
  drop_na()

# Where problems?
simu_id |>
  mutate(errors = map(simu_res, \(x) x$error)) |>
  unnest(errors) |> View()

# How many converged?
res |>
  summarise(
    count = sum(converged),
    .by = dist:method
  ) |>
  pivot_wider(names_from = method, values_from = count) |>
  print(n = Inf)

# Inspect twofac  n = 1000
res |>
  filter(model == "twofac", n == 1000, rel == 0.8) |>
  count(converged, dist, optim_message)
