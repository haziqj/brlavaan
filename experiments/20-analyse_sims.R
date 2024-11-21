library(tidyverse)
theme_set(theme_bw())
library(gt)
library(latex2exp)
load("experiments/simu_res_growth.RData")
load("experiments/simu_res_twofac.RData")

# Any errors?
# which(sapply(simu_res_growth, \(x) length(x$errors) != 0))
# x <- which(sapply(simu_res_twofac, \(x) length(x$errors) != 0))
# sapply(simu_res_twofac, \(x) nrow(x$simu_res))
# lapply(simu_res_twofac[x], \(x) length(x$errors))

res <-
  do.call(
    rbind,
    lapply(c(simu_res_growth, simu_res_twofac), \(x) x$simu_res)
  ) |>
  drop_na() |>
  mutate(
    model = factor(model, labels = c("Growth model", "Two factor model")),
    estimator = factor(estimator, levels = c("ML", "eBRM", "iBRM", "iBRMp")),
    dist = factor(dist, levels = c("Normal", "Kurtosis", "Non-normal")),
    rel = factor(rel, levels = c(0.8, 0.5), labels = c("Rel = 0.8", "Rel = 0.5")),
    n = factor(n),
    covered = truth <= est + qnorm(0.975) * se & truth >= est - qnorm(0.975) * se
  )

## ----- Table 1 ---------------------------------------------------------------
# Convergence failures Bootstrap resamples: Mean, median, min and max number of
# failed Bootstrap resamples for a single simulation (max = total number of
# resamples = 500).

res |>
  summarise(
    fail = any(!converged),
    .by = c(simu, dist, model, n, rel, estimator)
  ) |>
  summarise(
    count = sum(fail),
    .by = c(model, rel, n, estimator, dist)
  ) |>
  pivot_wider(names_from = c(dist, estimator), values_from = count) |>
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
    ends_with("eBRM") ~ "eBRM",
    ends_with("iBRM") ~ "iBRM",
    ends_with("iBRMp") ~ "iBRMp"
  ) |>
  tab_caption(md("`nlminb` non-convergence counts (maximum = 1000)"))
