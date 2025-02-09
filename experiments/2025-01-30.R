library(tidyverse)
library(brlavaan)
mycols <- c(
  ML = "#E31A1C",
  eRBM = "#A6CEE3",
  iRBM = "#1F78B4",
  Jackknife = "#B2DF8A",
  Bootstrap = "#33A02C",
  `Ozenne et al.` = "#FDBF6F",
  REML = "#FF7F00"
)

load("experiments/simu_res_growth.RData")
load("experiments/simu_res_twofac.RData")

# Create a simid table
simu_res <- c(simu_res_growth, simu_res_twofac)
simu_id <-
  expand_grid(
    model = c("growth", "twofac"),
    dist = c("Normal", "Kurtosis", "Non-normal"),
    rel = as.character(c(0.8, 0.5)),
    n = c(50, 100)
  ) |>
  rownames_to_column(var = "simid") |>
  mutate(simid = as.integer(simid))

# Create (nested) results data frame, used for counting convergences etc.
res_nested <-
  simu_res |>
  imap(\(x, idx) bind_cols(simid = idx, x$simu_res)) |>
  bind_rows() |>
  drop_na()

# This is the full (raw) results data frame
res <-
  res_nested |>
  mutate(param = lapply(truth, names), .before = est) |>
  unnest(param:truth) |>
  select(!starts_with("info")) |>
  # For the growth model, keep the first instance of param == "v"
  filter(
    row_number() == 1,
    .by = c(simid, sim, dist, model, rel, n, method, param)
  ) |>
  mutate(
    dist = factor(dist, levels = c("Normal", "Kurtosis", "Non-normal")),
    rel = factor(rel, levels = c("0.8", "0.5"), labels = c("Rel = 0.8", "Rel = 0.5")),
    method = factor(method, levels = rev(names(mycols))),
    n = as.factor(n)
  )

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
  select(-starts_with("info")) |>
  print(n = Inf)

res_nested |>
  count(model, method, optim_message)



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
  mutate(n = factor(n)) |>
  ggplot(aes(as.numeric(n), count, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(dist ~ model * rel) +
  scale_fill_manual(values = mycols) +
  scale_x_continuous(
    breaks = 1:2,
    labels = c(50, 100),
    name = "Sample size"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


twofacpars <- c("y1~~y1", "fx~~fx", "fy~~fy", "fy~fx", "fx=~x2")
growthpars <- c("v", "i~~i", "s~~s", "i~~s")

# fig_twofac_perf <-
res |>
  # mutate(converged = case_when(
  #   optim_message == "false convergence (8)" ~ TRUE,
  #   max_loglik > 0 ~ FALSE,
  #   TRUE ~ converged
  # )) |>
  # mutate(converged = all(converged), .by = c(simid, sim)) |>
  # filter(converged) |>
  mutate(
    bias = est - truth,
    relbias = bias / truth,
    covered = truth <= est + qnorm(0.975) * se & truth >= est - qnorm(0.975) * se
  ) |>
  filter(n == 50, dist == "Normal", param %in% twofacpars) |>
  summarise(
    mean_bias = mean(bias, na.rm = TRUE, trim = 0.05),
    median_bias = median(bias, na.rm = TRUE),
    rmse = sqrt(mean(bias ^ 2, na.rm = TRUE, trim = 0.05)),
    truth = first(truth),
    coverage = mean(covered, na.rm = TRUE),
    .by = c(param, rel, method, dist, n)
  ) |>
  mutate(
    mean_bias = mean_bias / truth,
    rmse = rmse / truth
  ) |>
  pivot_longer(
    cols = c(mean_bias, rmse, coverage),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(metric = factor(
    metric,
    levels = c("mean_bias",  "rmse", "coverage"),
    labels = c("Relative mean bias",  "Relative RMSE", "Coverage")
  )) |>
  ggplot(aes(value, fct_rev(param), fill = fct_rev(method))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_vline(
    data = tibble(
      metric = factor(c("Relative mean bias", "Relative RMSE", "Coverage")),
      value = c(0, 0, 0.95),
    ),
    aes(xintercept = value),
    linetype = "dashed"
  ) +
  facet_grid(rel * dist ~ metric, scales = "free_x") +
  ggh4x::facetted_pos_scales(
    x = list(
      scale_x_continuous(labels = scales::percent),
      scale_x_continuous(labels = scales::percent),
      scale_x_continuous(limits = c(0.55, 1), labels = scales::percent)
    )
  ) +
  # scale_y_discrete(labels = rev(c(
  #   expression(theta),
  #   expression(Psi["1,1"]),
  #   expression(Psi["2,2"]),
  #   expression(beta)
  #   # expression(Lambda["2,1"])
  # ))) +
  scale_fill_manual(values = rev(mycols))  +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(
    x = NULL,
    y = NULL,
    fill = NULL,
    title = "n = 50"
  ) +
  theme_bw()


