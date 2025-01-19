# source("experiments/20-analyse_sims.R")

# From D&R
load("experiments/2FSEM_est_combined_final.RData")
res_dr <-
  as_tibble(Results) |>
  unite("simu", jobid, iteration, sep = ".") |>
  mutate(simu = as.numeric(factor(simu))) |>
  select(
    simu,
    method,
    dist,
    n = nobs,
    rel,
    `est_y1~~y1`,
    `est_fx~~fx`,
    `est_fy~~fy`,
    `est_fy~fx`,
    `est_fx=~x2`,
    `se_y1~~y1`,
    `se_fx~~fx`,
    `se_fy~~fy`,
    `se_fy~fx`,
    `se_fx=~x2`
  ) |>
  pivot_longer(
    cols = c(starts_with("est_"), starts_with("se_")),
    names_to = c(".value", "param"),
    names_sep = "_"
  ) |>
  filter(method %in% c("JB", "BB", "Ozenne", "REML")) |>
  # Some duplicates in the BB method??
  distinct(simu, method, dist, n, rel, param, est, se)

levels(res_dr$dist) <- c("Kurtosis", "Non-normal", "Normal")
res_dr$dist <- factor(res_dr$dist, levels = c("Normal", "Kurtosis", "Non-normal"))
levels(res_dr$rel) <- c("Rel = 0.5", "Rel = 0.8")
res_dr$rel <- factor(res_dr$rel, levels = c("Rel = 0.8", "Rel = 0.5"))

tmp_lev <- levels(res_dr$method)
tmp_lev <- gsub("JB", "Jackknife", tmp_lev)
tmp_lev <- gsub("BB", "Bootstrap", tmp_lev)
tmp_lev <- gsub("Ozenne", "Ozenne et al.", tmp_lev)
tmp_lev <- gsub("REML", "REML", tmp_lev)
levels(res_dr$method) <- tmp_lev

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
  select(simu, param, method, dist, n, rel, est, truth, covered) |>
  bind_rows(res_dr) |>
  mutate(param = factor(param, levels = twofacpars)) |>
  summarise(
    mean_bias = mean(est - truth, na.rm = TRUE, trim = 0.05),
    median_bias = median(est - truth, na.rm = TRUE),
    rmse = sqrt(mean((est - truth) ^ 2, na.rm = TRUE, trim = 0.05)),
    truth = first(truth),
    coverage = mean(covered, na.rm = TRUE, trim = 0.05),
    .by = c(param, rel, method, dist, n)
  ) |>
  mutate(
    rel_mean_bias = mean_bias / truth,
    rel_median_bias = median_bias / truth
  ) |>
  select(param:n, rel_mean_bias, rel_median_bias, rmse, coverage) |>
  arrange(param, method)

tab5 <-
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
    "method" ~ "Estimator"
  ) |>
  fmt_number(decimals = 2) |>
  fmt_markdown(columns = 1) |>
  cols_align(align = "center") |>
  tab_caption("Reliability = 0.8")

## ----- Table 6 ---------------------------------------------------------------
tab6 <-
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
    "method" ~ "Estimator"
  ) |>
  fmt_number(decimals = 2) |>
  fmt_markdown(columns = 1) |>
  cols_align(align = "center") |>
  tab_caption("Reliability = 0.5")

## ----- Figure 6 --------------------------------------------------------------
fig6 <-
  tab_twofac |>
  filter(dist != "Kurtosis") |>
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
  ggplot(aes(x = as.numeric(n), y = rel_median_bias, col = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "gray30") +
  geom_line() +
  geom_point(size = 1) +
  facet_grid(param ~ rel * dist, labeller = labeller(param = label_parsed)) +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000)) +
  scale_y_continuous(labels = scales::percent) +
  scale_colour_manual(values = mycols) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(
    x = "Sample size (n)",
    y = "Relative median bias",
    col = NULL
  ) +
  theme(legend.position = "top")

## ----- Figure 7 --------------------------------------------------------------
fig7 <-
  tab_twofac |>
  filter(dist != "Kurtosis") |>
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
  ggplot(aes(x = as.numeric(n), y = rel_mean_bias, col = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "gray30") +
  geom_line() +
  geom_point(size = 1) +
  facet_grid(param ~ rel * dist, labeller = labeller(param = label_parsed)) +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000)) +
  scale_y_continuous(labels = scales::percent) +
  scale_colour_manual(values = mycols) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(
    x = "Sample size (n)",
    y = "Relative mean bias",
    col = NULL
  ) +
  # coord_cartesian(ylim = c(-0.4, 0.3)) +
  theme(legend.position = "top")

## ----- Figure 8 --------------------------------------------------------------
fig8 <-
  left_join(
    tab_twofac,
    tab_twofac |>
      filter(method == "ML") |>
      select(param, rel, dist, n, rmse_ML = rmse)
  ) |>
  mutate(
    rel_rmse = rmse / rmse_ML,
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
  # filter(dist != "Kurtosis") |>
  ggplot(aes(x = as.numeric(n), y = rel_rmse, col = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "gray30") +
  geom_line() +
  geom_point(size = 1) +
  facet_grid(param ~ rel * dist, labeller = labeller(param = label_parsed)) +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000)) +
scale_colour_manual(values = mycols) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(
    x = "Sample size (n)",
    y = "Relative RMSE (with respect to ML)",
    col = NULL
  ) +
  theme(legend.position = "top")

## ----- Table 7 ---------------------------------------------------------------
tab7 <-
  tab_twofac |>
  mutate(param = factor(param, labels = c(
    "$\\theta$",
    "$\\Psi_{1,1}$",
    "$\\Psi_{2,2}$",
    "$\\beta$",
    "$\\Lambda_{2,1}$"
  ))) |>
  filter(dist != "Kurtosis") |>
  pivot_wider(
    id_cols = c(param, method),
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
    "method" ~ "Estimator"
  ) |>
  fmt_number(decimals = 2) |>
  fmt_markdown(columns = 1) |>
  cols_align(align = "center")

## ----- New plots -------------------------------------------------------------
ress <-
  res |>
  filter(model == "Two factor model", param %in% twofacpars) |>
  select(simu, param, method, dist, n, rel, est, truth, covered) |>
  bind_rows(res_dr) |>
  mutate(
    param = factor(param, levels = twofacpars),
    bias = est - truth,
    relbias = bias / truth,
    n = as.numeric(as.character(n)),
    nbias = n * bias,
    highlight = case_when(
      grepl("RBM", method) ~ "RBM",
      TRUE ~ "Non-RBM"
    ),
    under = est < truth
  )

# Distribution plot
fig_twofac_dis <-
  ress |>
  mutate(method = fct_rev(method)) |>
  filter(n == 50, dist != "Kurtosis", param != "fx=~x2") |>
  ggplot(aes(relbias, fct_rev(param), fill = method)) +
  geom_point(alpha = 0.1, size = 0.3, aes(group = method),
             position = position_dodge(width = .75)) +
  geom_boxplot(outliers = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = rev(mycols)) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_discrete(labels = rev(c(
    expression(theta),
    expression(Psi["1,1"]),
    expression(Psi["2,2"]),
    expression(beta)
    # expression(Lambda["2,1"])
  ))) +
  coord_cartesian(xlim = c(-3, 3)) +
  guides(fill = guide_legend(reverse = TRUE, nrow = 1)) +
  facet_grid(rel ~ dist, scales = "free_y") +
  labs(
    x = "Relative bias",
    y = NULL,
    fill = NULL
  ) +
  theme_bw() +
  theme(legend.position = "top"); fig_twofac_dis

# Performance plot
fig_twofac_perf <-
  ress |>
  filter(n == 50, dist != "Kurtosis", param != "fx=~x2") |>
  summarise(
    mean_bias = mean(est - truth, na.rm = TRUE, trim = 0.05),
    median_bias = median(est - truth, na.rm = TRUE),
    rmse = sqrt(mean((est - truth) ^ 2, na.rm = TRUE, trim = 0.05)),
    truth = first(truth),
    coverage = mean(covered, na.rm = TRUE),
    pu = mean(under, na.rm = TRUE),
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
  scale_y_discrete(labels = rev(c(
    expression(theta),
    expression(Psi["1,1"]),
    expression(Psi["2,2"]),
    expression(beta)
    # expression(Lambda["2,1"])
  ))) +
  scale_fill_manual(values = rev(mycols))  +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_bw(); fig_twofac_perf


