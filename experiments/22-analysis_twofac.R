# source("experiments/20-analyse_sims.R")
# From D&R
load("experiments/2FSEM_est_combined_final.RData")
x <-
  as_tibble(Results) |>
  select(
    simu = iteration,
    estimator = method,
    dist,
    n = nobs,
    rel,
    `est_y1~~y1`,
    `est_fx~~fx`,
    `est_fy~~fy`,
    `est_fy~fx`,
    `est_fx=~x2`
  ) |>
  mutate(id = row_number()) |>
  pivot_longer(
    cols = starts_with("est_"),
    names_to = "param",
    values_to = "est"
  ) |>
  mutate(param = gsub("est_", "", param)) |>
  filter(estimator %in% c("JB", "BB", "Ozenne", "REML"))
y <-
  as_tibble(Results) |>
  select(
    simu = iteration,
    estimator = method,
    dist,
    n = nobs,
    rel,
    `se_y1~~y1`,
    `se_fx~~fx`,
    `se_fy~~fy`,
    `se_fy~fx`,
    `se_fx=~x2`
  ) |>
  mutate(id = row_number()) |>
  pivot_longer(
    cols = starts_with("se_"),
    names_to = "param",
    values_to = "se"
  ) |>
  mutate(param = gsub("se_", "", param)) |>
  filter(estimator %in% c("JB", "BB", "Ozenne", "REML"))
res_dr <- left_join(x, y) |> select(-id)
rm(x)
rm(y)

levels(res_dr$dist) <- c("Kurtosis", "Non-normal", "Normal")
res_dr$dist <- factor(res_dr$dist, levels = c("Normal", "Kurtosis", "Non-normal"))
res_dr$param <- gsub("est_", "", res_dr$param)
levels(res_dr$rel) <- c("Rel = 0.5", "Rel = 0.8")
res_dr$rel <- factor(res_dr$rel, levels = c("Rel = 0.8", "Rel = 0.5"))
levels(res_dr$estimator)[c(2, 3, 4, 8)] <- c(
  "Jackknife", "Bootstrap", "Ozenne et al.", "REML"
)

res_dr <-
  left_join(
    res_dr,
    summarise(
      filter(res, model == "Two factor model"),
      truth = first(truth),
      .by = c(param, rel)
    )
  ) |>
  mutate(
    covered = truth <= est + qnorm(0.975) * se & truth >= est - qnorm(0.975) * se
  ) |>
  select(-se)


## ----- Table 5 ---------------------------------------------------------------
# Results (trimmed) two-factor model reliability .80.

twofacpars <- c("y1~~y1", "fx~~fx", "fy~~fy", "fy~fx", "fx=~x2")

tab_twofac <-
  res |>
  filter(model == "Two factor model") |>
  filter(param %in% twofacpars) |>
  select(simu, param, estimator, dist, n, rel, est, truth, covered) |>
  # bind_rows(res_dr) |>
  mutate(param = factor(param, levels = twofacpars)) |>
  summarise(
    mean_bias = mean(est - truth, na.rm = TRUE, trim = 0.05),
    median_bias = median(est - truth, na.rm = TRUE),
    rmse = sqrt(mean((est - truth) ^ 2, na.rm = TRUE, trim = 0.05)),
    truth = first(truth),
    coverage = mean(covered, na.rm = TRUE, trim = 0.05),
    .by = c(param, rel, estimator, dist, n)
  ) |>
  mutate(
    rel_mean_bias = mean_bias / truth,
    rel_median_bias = median_bias / truth
  ) |>
  select(param:n, rel_mean_bias, rel_median_bias, rmse, coverage) |>
  arrange(param, estimator)

# gt table
tab_twofac |>
  mutate(param = factor(param, labels = c(
    "$\\theta$",
    "$\\Psi_{1,1}$",
    "$\\Psi_{2,2}$",
    "$\\beta$",
    "$\\Lambda_{2,1}$"
  ))) |>
  filter(dist != "Kurtosis") |>
  select(-coverage) |>
  pivot_wider(
    names_from = c(dist, n),
    values_from = c(rel_mean_bias, rel_median_bias, rmse)
  ) |>
  filter(rel == "Rel = 0.8") |>  # RELIABILITY
  select(-rel) |>
  gt(
    groupname_col = c("param"),
    process_md = TRUE
  ) |>
  tab_spanner(
    label = "Relative mean bias",
    columns = contains("rel_mean_bias")
  ) |>
  tab_spanner(
    label = "Relative median bias",
    columns = contains("rel_median_bias")
  ) |>
  tab_spanner(
    label = "RMSE",
    columns = contains("rmse")
  ) |>
  tab_spanner(
    label = "Normal data",
    columns = contains("Normal", ignore.case = FALSE)
  ) |>
  tab_spanner(
    label = "Non-normal data",
    columns = contains("Non-normal")
  ) |>
  cols_label(
    ends_with("15") ~ "15",
    ends_with("20") ~ "20",
    ends_with("50") ~ "50",
    ends_with("100") ~ "100",
    ends_with("1000") ~ "1000",
    "estimator" ~ "Estimator"
  ) |>
  fmt_number(decimals = 2) |>
  fmt_markdown(columns = 1) |>
  cols_align(align = "center") |>
  tab_caption("Reliability = 0.8")

## ----- Table 6 ---------------------------------------------------------------
tab_twofac |>
  mutate(param = factor(param, labels = c(
    "$\\theta$",
    "$\\Psi_{1,1}$",
    "$\\Psi_{2,2}$",
    "$\\beta$",
    "$\\Lambda_{2,1}$"
  ))) |>
  filter(dist != "Kurtosis") |>
  select(-coverage) |>
  pivot_wider(
    names_from = c(dist, n),
    values_from = c(rel_mean_bias, rel_median_bias, rmse)
  ) |>
  filter(rel == "Rel = 0.5") |>
  select(-rel) |>
  gt(
    groupname_col = c("param"),
    process_md = TRUE
  ) |>
  tab_spanner(
    label = "Relative mean bias",
    columns = contains("rel_mean_bias")
  ) |>
  tab_spanner(
    label = "Relative median bias",
    columns = contains("rel_median_bias")
  ) |>
  tab_spanner(
    label = "RMSE",
    columns = contains("rmse")
  ) |>
  tab_spanner(
    label = "Normal data",
    columns = contains("Normal", ignore.case = FALSE)
  ) |>
  tab_spanner(
    label = "Non-normal data",
    columns = contains("Non-normal")
  ) |>
  cols_label(
    ends_with("15") ~ "15",
    ends_with("20") ~ "20",
    ends_with("50") ~ "50",
    ends_with("100") ~ "100",
    ends_with("1000") ~ "1000",
    "estimator" ~ "Estimator"
  ) |>
  fmt_number(decimals = 2) |>
  fmt_markdown(columns = 1) |>
  cols_align(align = "center") |>
  tab_caption("Reliability = 0.5")

## ----- Figure 3 --------------------------------------------------------------
tab_twofac |>
  mutate(
    param = factor(
      param,
      labels = c(
        expression(theta),
        expression(Psi["1,1"]),
        expression(Psi["2,2"]),
        expression(beta),
        expression(Lambda["2,1"])
      )
    )
  ) |>
  filter(dist != "Kurtosis") |>
  ggplot(aes(x = as.numeric(n), y = rel_median_bias, col = estimator)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "gray30") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1) +
  facet_grid(param ~ rel * dist, labeller = labeller(param = label_parsed)) +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000)) +
  ggsci::scale_colour_d3() +
  labs(
    x = "Sample size (n)",
    y = "Relative median bias",
    col = NULL
  ) +
  # coord_cartesian(ylim = c(-0.4, 0.3)) +
  theme(legend.position = "top")

## ----- Figure 4 --------------------------------------------------------------
tab_growth |>
  mutate(
    param = factor(
      param,
      labels = c(
        expression(theta["1,1"]),
        expression(Psi["1,1"]),
        expression(Psi["2,2"]),
        expression(Psi["1,2"])
      )
    )
  ) |>
  filter(dist != "Kurtosis") |>
  ggplot(aes(x = as.numeric(n), y = rel_mean_bias, col = estimator)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "gray30") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1) +
  facet_grid(param ~ rel * dist, labeller = labeller(param = label_parsed)) +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000)) +
  ggsci::scale_colour_d3() +
  labs(
    x = "Sample size (n)",
    y = "Relative mean bias",
    col = NULL
  ) +
  coord_cartesian(ylim = c(-0.4, 0.3)) +
  theme(legend.position = "top")

## ----- Figure 5 --------------------------------------------------------------
left_join(
  tab_growth,
  tab_growth |>
    filter(estimator == "ML") |>
    select(param, rel, dist, n, rmse_ML = rmse)
) |>
  mutate(
    rel_rmse = rmse / rmse_ML,
    param = factor(
      param,
      labels = c(
        expression(theta["1,1"]),
        expression(Psi["1,1"]),
        expression(Psi["2,2"]),
        expression(Psi["1,2"])
      )
    )
  ) |>
  filter(dist != "Kurtosis") |>
  ggplot(aes(x = as.numeric(n), y = rel_rmse, col = estimator)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "gray30") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1) +
  facet_grid(param ~ rel * dist, labeller = labeller(param = label_parsed)) +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000)) +
  ggsci::scale_colour_d3() +
  labs(
    x = "Sample size (n)",
    y = "Relative RMSE (with respect to ML)",
    col = NULL
  ) +
  coord_cartesian(ylim = c(0, 1.1)) +
  theme(legend.position = "top")


## table
tab_growth |>
  mutate(
    param = factor(param, labels = c("$\\theta_{1,1}$", "$\\Psi_{1,1}$", "$\\Psi_{2,2}$", "$\\Psi_{1,2}$"))
  ) |>
  filter(dist != "Kurtosis") |>
  pivot_wider(
    id_cols = c(param, estimator),
    names_from = c(dist, rel, n),
    values_from = c(coverage)
  ) |>
  gt(
    groupname_col = c("param"),
    process_md = TRUE
  ) |>
  tab_spanner(
    label = "Reliability = 0.80",
    columns = contains("Rel = 0.8")
  ) |>
  tab_spanner(
    label = "Reliability = 0.50",
    columns = contains("Rel = 0.5")
  ) |>
  tab_spanner(
    label = "Normal data",
    columns = contains("Normal", ignore.case = FALSE)
  ) |>
  tab_spanner(
    label = "Non-normal data",
    columns = contains("Non-normal")
  ) |>
  cols_label(
    ends_with("15") ~ "15",
    ends_with("20") ~ "20",
    ends_with("50") ~ "50",
    ends_with("100") ~ "100",
    ends_with("1000") ~ "1000",
    "estimator" ~ "Estimator"
  ) |>
  fmt_number(decimals = 2) |>
  fmt_markdown(columns = 1) |>
  cols_align(align = "center")


















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
  geom_point(size = 1) +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000)) +
  ggsci::scale_colour_d3() +
  facet_grid(param ~ rel + dist) +
  labs(
    x = "Sample size (n)",
    y = "Relative median bias",
    col = NULL
  ) +
  theme(legend.position = "top"); p_twofac

save(p_growth, p_twofac, file = "R/DR_sims.RData")
