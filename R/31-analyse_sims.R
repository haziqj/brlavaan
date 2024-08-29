library(tidyverse)
theme_set(theme_bw())
load("R/sim_results.RData")

# Overall plot
results |>
  # drop_na() |>
  filter(abs(est - truth) / truth < 1) |>
  # filter(method %in% c("iRBM", "ML")) |>
  group_by(dist, model, rel, n, method) |>
  summarise(bias = mean(est - truth, na.rm = TRUE)) |>
  ggplot(aes(x = as.numeric(n), y = abs(bias), col = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "gray30") +
  geom_line(size = 0.8) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", \(x) 10 ^ x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000)) +
  facet_grid(model + rel ~ dist, scales = "free_y") +
  labs(
    x = "Sample size (n)",
    y = "Absolute mean bias",
    col = NULL
  ) +
  theme(legend.position = "top")

# Two factor model results -----------------------------------------------------
results |>
  filter(model == "Two factor model", dist %in% c("Normal", "Non-normal")) |>
  filter(param %in% c("y1~~y1", "fx~~fx", "fy~~fy", "fy~fx", "fx=~x2")) |>
  filter(abs(est) < 2) |>
  group_by(param, dist, rel, n, method) |>
  summarise(
    bias = mean(est - truth, na.rm = TRUE),
    truth = first(truth)
  ) |>
  mutate(
    rel_bias = bias / truth,
    param = factor(param, levels = c("y1~~y1", "fx~~fx", "fy~~fy", "fy~fx", "fx=~x2")),
  ) |>
  ggplot(aes(x = as.numeric(n), y = rel_bias, col = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "gray30") +
  geom_line() +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000)) +
  facet_grid(param ~ rel + dist, scales = "free_y") +
  labs(
    x = "Sample size (n)",
    y = "Relative mean bias",
    col = NULL
  ) +
  theme(legend.position = "top")

results |>
  filter(model == "Two factor model") |>
  filter(param %in% c("y1~~y1", "fx~~fx", "fy~~fy", "fy~fx", "fx=~x2")) |>
  filter(abs(est) < 5) |>
  group_by(dist, rel, param, n, method) |>
  summarise(
    # bias = mean(est - truth, na.rm = TRUE),
    mse = mean((est - truth) ^ 2, na.rm = TRUE),
    # cover = mean((truth - 1.96 * se < est) & (est < truth + 1.96 * se), na.rm = TRUE),
    # time = mean(time, na.rm = TRUE),
    # nsimu = n()
  ) |>
  pivot_wider(names_from = n, values_from = mse, names_prefix = "n = ") |>
  print(n = 100)

