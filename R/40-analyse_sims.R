library(tidyverse)
theme_set(theme_bw())
load("R/sim_results.RData")
load("R/sim_results_twofac.RData")

# Growth results ---------------------------------------------------------------
p_growth <-
  results |>
  filter(model == "Growth model", dist %in% c("Normal", "Non-normal")) |>
  filter(param %in% c("v", "i~~i", "s~~s", "i~~s")) |>
  # filter(abs(est - truth) / truth < 1) |>
  filter(method != "iRBMp") |>
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
  geom_point(size = 1) +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000)) +
  facet_grid(param ~ rel + dist) +
  labs(
    x = "Sample size (n)",
    y = "Relative median bias",
    col = NULL
  ) +
  theme(legend.position = "top"); p_growth

# Two factor model results -----------------------------------------------------
p_twofac <-
  res_twofac |>
  filter(dist %in% c("Normal", "Non-normal")) |>
  filter(param %in% c("y1~~y1", "fx~~fx", "fy~~fy", "fy~fx", "fx=~x2")) |>
  # filter(method != "iRBMp") |>
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
  facet_grid(param ~ rel + dist) +
  labs(
    x = "Sample size (n)",
    y = "Relative median bias",
    col = NULL
  ) +
  theme(legend.position = "top"); p_twofac

save(p_growth, p_twofac, file = "R/DR_sims.RData")
