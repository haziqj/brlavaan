# source("experiments/20-analyse_sims.R")

## ----- Table 2 ---------------------------------------------------------------
# Results (trimmed) growth curve model reliability .80.

res_v <-
  res |>
  filter(model == "Growth model") |>
  filter(param %in% c("v")) |>
  summarise(
    across(c(est, truth, covered), first),
    .by = c(simu, param, method, dist, n, rel)
  )
res_nonv <- res |>
  filter(model == "Growth model") |>
  filter(param %in% c("i~~i", "s~~s", "i~~s")) |>
  select(simu, param, method, dist, n, rel, est, truth, covered)

truth <-
  bind_rows(
    res_v,
    res_nonv
  ) |>
  summarise(
    truth = first(truth),
    .by = c(param, rel)
  )

# From D&R
load("experiments/GCM_est_combined_final.RData")
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
    `est_Day0~~Day0`,
    `est_i~~i`,
    `est_s~~s`,
    `est_i~~s`,
    `se_Day0~~Day0`,
    `se_i~~i`,
    `se_s~~s`,
    `se_i~~s`
  ) |>
  pivot_longer(
    cols = c(starts_with("est_"), starts_with("se_")),
    names_to = c(".value", "param"),
    names_sep = "_"
  ) |>
  filter(method %in% c("JB", "BB", "Ozenne", "REML")) |> # NO NEED REML
  distinct(simu, method, dist, n, rel, param, est, se)


levels(res_dr$dist) <- c("Kurtosis", "Non-normal", "Normal")
res_dr$dist <- factor(res_dr$dist, levels = c("Normal", "Kurtosis", "Non-normal"))
res_dr$param <- gsub("est_", "", res_dr$param)
res_dr$param <- gsub("Day0~~Day0", "v", res_dr$param)
levels(res_dr$rel) <- c("Rel = 0.5", "Rel = 0.8")
res_dr$rel <- factor(res_dr$rel, levels = c("Rel = 0.8", "Rel = 0.5"))

tmp_lev <- levels(res_dr$method)
tmp_lev <- gsub("JB", "Jackknife", tmp_lev)
tmp_lev <- gsub("BB", "Bootstrap", tmp_lev)
tmp_lev <- gsub("Ozenne", "Ozenne et al.", tmp_lev)
tmp_lev <- gsub("REML", "REML", tmp_lev)
levels(res_dr$method) <- tmp_lev

res_dr <-
  left_join(res_dr, truth) |>
  mutate(
    covered = truth <= est + qnorm(0.975) * se & truth >= est - qnorm(0.975) * se
  ) |>
  select(-se)

tab_growth <-
  bind_rows(
    res_v,
    res_nonv,
    res_dr
  ) |>
  mutate(param = factor(param, levels = c("v", "i~~i", "s~~s", "i~~s"))) |>
  summarise(
    mean_bias = mean(est - truth, na.rm = TRUE, trim = 0.05),
    median_bias = median(est - truth, na.rm = TRUE),
    rmse = sqrt(mean((est - truth) ^ 2, na.rm = TRUE, trim = 0.05)),
    truth = first(truth),
    coverage = mean(covered, na.rm = TRUE),
    .by = c(param, rel, method, dist, n)
  ) |>
  mutate(
    rel_mean_bias = mean_bias / truth,
    rel_median_bias = median_bias / truth
  ) |>
  select(param:n, rel_mean_bias, rel_median_bias, rmse, coverage) |>
  arrange(param, method)

tab2 <-
  tab_growth |>
  mutate(
    param = factor(param, labels = c("$\\theta_{1,1}$", "$\\Psi_{1,1}$", "$\\Psi_{2,2}$", "$\\Psi_{1,2}$"))
  ) |>
  filter(dist != "Kurtosis") |>
  select(-coverage) |>
  pivot_wider(
    names_from = c(dist, n),
    values_from = c(rel_mean_bias, rel_median_bias, rmse)
  ) |>
  filter(rel == "Rel = 0.8") |>
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
    "method" ~ "method"
  ) |>
  fmt_number(decimals = 2) |>
  fmt_markdown(columns = 1) |>
  cols_align(align = "center")

## ----- Table 3 ---------------------------------------------------------------
tab3 <-
  tab_growth |>
  mutate(
    param = factor(param, labels = c("$\\theta_{1,1}$", "$\\Psi_{1,1}$", "$\\Psi_{2,2}$", "$\\Psi_{1,2}$"))
  ) |>
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
    "method" ~ "method"
  ) |>
  fmt_number(decimals = 2) |>
  fmt_markdown(columns = 1) |>
  cols_align(align = "center")

## ----- Figure 3 --------------------------------------------------------------
fig3 <-
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
  ggplot(aes(x = as.numeric(n), y = rel_median_bias, col = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "gray30") +
  geom_point(size = 1) +
  geom_line() +
  facet_grid(param ~ rel * dist, labeller = labeller(param = label_parsed)) +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000)) +
  scale_y_continuous(labels = scales::percent) +
  scale_colour_manual(values = mycols) +
  labs(
    x = "Sample size (n)",
    y = "Relative median bias",
    col = NULL
  ) +
  coord_cartesian(ylim = c(-0.4, 0.3)) +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.position = "top")

## ----- Figure 4 --------------------------------------------------------------
fig4 <-
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
  # filter(method %in% c("ML", "eRBM", "iRBM")) |>
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
  coord_cartesian(ylim = c(-0.4, 0.3)) +
  theme(legend.position = "top"); fig4

## ----- Figure 5 --------------------------------------------------------------
fig5 <-
  left_join(
    tab_growth,
    tab_growth |>
      filter(method == "ML") |>
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
  coord_cartesian(ylim = c(0.9, 1.3)) +
  theme(legend.position = "top")

## ----- Table 4 ---------------------------------------------------------------
tab4 <-
  tab_growth |>
  mutate(
    param = factor(param, labels = c("$\\theta_{1,1}$", "$\\Psi_{1,1}$", "$\\Psi_{2,2}$", "$\\Psi_{1,2}$"))
  ) |>
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
    "method" ~ "method"
  ) |>
  fmt_number(decimals = 2) |>
  fmt_markdown(columns = 1) |>
  cols_align(align = "center")

## ----- New plots -------------------------------------------------------------
ress <-
  bind_rows(
    res_v,
    res_nonv,
    res_dr
  ) |>
  mutate(
    param = factor(
      param,
      levels = c("v", "i~~i", "s~~s", "i~~s")
    ),
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
fig_growth_dis <-
  ress |>
  mutate(method = fct_rev(method)) |>
  filter(n == 50, dist != "Kurtosis", param != "v") |>
  ggplot(aes(relbias, fct_rev(param), fill = method)) +
  geom_point(alpha = 0.1, size = 0.3, aes(group = method),
             position = position_dodge(width = .75)) +
  geom_boxplot(outliers = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = rev(mycols)) +
  scale_x_continuous(labels = scales::percent) +
  scale_y_discrete(labels = rev(c(
    expression(theta["1,1"]),
    expression(Psi["1,1"]),
    expression(Psi["2,2"]),
    expression(Psi["1,2"])
  ))) +
  coord_cartesian(xlim = c(-2, 2)) +
  guides(fill = guide_legend(reverse = TRUE, nrow = 1)) +
  facet_grid(rel ~ dist, scales = "free_y") +
  labs(
    x = "Relative bias",
    y = NULL,
    fill = NULL
  ) +
  theme_bw() +
  theme(legend.position = "top")

# Performance plot
fig_growth_perf <-
  ress |>
  filter(n == 50, dist != "Kurtosis", param != "v") |>
  summarise(
    mean_bias = mean(est - truth, na.rm = TRUE, trim = 0.05),
    median_bias = median(est - truth, na.rm = TRUE),
    rmse = sqrt(mean((est - truth) ^ 2, na.rm = TRUE)),
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
    expression(theta["1,1"]),
    expression(Psi["1,1"]),
    expression(Psi["2,2"]),
    expression(Psi["1,2"])
  ))) +
  scale_fill_manual(values = rev(mycols))  +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(
    x = NULL,
    y = NULL,
    fill = NULL
  ) +
  theme_bw()


