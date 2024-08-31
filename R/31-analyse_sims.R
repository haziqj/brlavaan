library(tidyverse)
theme_set(theme_bw())
load("R/sim_results.RData")

# Growth results ---------------------------------------------------------------
p_growth <-
  results |>
  filter(model == "Growth model", dist %in% c("Normal", "Non-normal")) |>
  filter(param %in% c("v", "i~~i", "s~~s", "i~~s", "fx=~x2")) |>
  # filter(abs(est) < 1) |>
  group_by(param, dist, rel, n, method) |>
  summarise(
    bias = median(est - truth, na.rm = TRUE),
    truth = first(truth)
  ) |>
  mutate(
    rel_bias = bias / truth,
    param = factor(param, levels = c("v", "i~~i", "s~~s", "i~~s", "fx=~x2")),
  ) |>
  ggplot(aes(x = as.numeric(n), y = rel_bias, col = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "gray30") +
  geom_line(linewidth = 0.8) +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000)) +
  facet_grid(param ~ rel + dist, scales = "free_y") +
  labs(
    x = "Sample size (n)",
    y = "Relative median bias",
    col = NULL
  ) +
  theme(legend.position = "top")


# Two factor model results -----------------------------------------------------
p_twofac <-
  results |>
  filter(model == "Two factor model", dist %in% c("Normal", "Non-normal")) |>
  filter(param %in% c("y1~~y1", "fx~~fx", "fy~~fy", "fy~fx", "fx=~x2")) |>
  # filter(abs(est) < 1) |>
  group_by(param, dist, rel, n, method) |>
  summarise(
    bias = median(est - truth, na.rm = TRUE),
    truth = first(truth)
  ) |>
  mutate(
    rel_bias = bias / truth,
    param = factor(param, levels = c("y1~~y1", "fx~~fx", "fy~~fy", "fy~fx", "fx=~x2")),
  ) |>
  ggplot(aes(x = as.numeric(n), y = rel_bias, col = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "gray30") +
  geom_line(linewidth = 0.8) +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000)) +
  facet_grid(param ~ rel + dist, scales = "free_y") +
  labs(
    x = "Sample size (n)",
    y = "Relative median bias",
    col = NULL
  ) +
  theme(legend.position = "top")


save(p_growth, p_twofac, file = "R/DR_sims.RData")
