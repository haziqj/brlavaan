# Where problems?
simu_id |>
  mutate(
    errors = map(simu_res, \(x) x$error),
    errors = sapply(errors, as.character)
  ) |>
  filter(sapply(errors, length) > 0) |>
  unnest(errors)

# How many converged?
res_nested |>
  summarise(
    count = sum(converged),
    .by = dist:method
  ) |>
  pivot_wider(names_from = method, values_from = count) |>
  print(n = Inf)

res_nested |>
  # Here's where the non-convergences are filtered out
  # mutate(converged = case_when(
  #   optim_message == "false convergence (8)" ~ TRUE,
  #   TRUE ~ converged
  # )) |>
  # mutate(converged = all(converged), .by = c(simid, sim)) |>
  # filter(converged) |>
  summarise(
    count = sum(converged),
    .by = dist:method
  ) |>
  ggplot(aes(as.numeric(n), count, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(model ~ dist * rel) +
  scale_fill_manual(values = mycols) +
  scale_x_continuous(breaks = 1:5, labels = c(15, 20, 50, 100, 1000),
                     name = "Sample size")

## ----- Table 1 ---------------------------------------------------------------
# Convergence failures Bootstrap resamples: Mean, median, min and max number of
# failed Bootstrap resamples for a single simulation (max = total number of
# resamples = 500).

tab1 <-
  res |>
  summarise(
    fail = any(!converged),
    .by = c(simu, dist, model, n, rel, method)
  ) |>
  summarise(
    count = sum(!fail),
    .by = c(model, rel, n, method, dist)
  ) |>
  pivot_wider(names_from = c(dist, method), values_from = count) |>
  gt(
    rowname_col = "n",
    groupname_col = c("model", "rel")
  ) |>
  tab_spanner(
    label = "Normal",
    columns = starts_with("Normal")
  ) |>
  tab_spanner(
    label = "Kurtosis",
    columns = starts_with("Kurtosis")
  ) |>
  tab_spanner(
    label = "Non-normal",
    columns = starts_with("Non-normal")
  ) |>
  cols_label(
    ends_with("ML") ~ "ML",
    ends_with("eRBM") ~ "eBR",
    ends_with("iRBM") ~ "iBR"
  ) |>
  tab_caption(md("`nlminb` non-convergence counts (maximum = 1000)"))

timing <-
  res |>
  summarise(
    timing = mean(timing),
    .by = c(model, n, method)
  ) |>
  pivot_wider(names_from = method, values_from = timing)
