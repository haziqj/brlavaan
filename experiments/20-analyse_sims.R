library(tidyverse)
theme_set(theme_bw())
library(gt)
library(latex2exp)
load(here::here("experiments/simu_res_growth.RData"))
load(here::here("experiments/simu_res_twofac.RData"))

## ----- Download D&R sims -----------------------------------------------------
dr_file1 <- here::here("experiments/GCM_est_combined_final.RData")
dr_file2 <- here::here("experiments/2FSEM_est_combined_final.RData")
if (!file.exists(dr_file1))
  download.file("https://osf.io/vjq5m/download", destfile = dr_file1)
if (!file.exists(dr_file2))
  download.file("https://osf.io/cw5b7/download", destfile = dr_file2)

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

tab1 <-
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

## ----- Growth & Two-factor results -------------------------------------------
source(here::here("experiments/21-analysis_growth.R"))
source(here::here("experiments/22-analysis_twofac.R"))

## ----- Save ------------------------------------------------------------------
save(
  tab1, tab2, tab3, tab4, tab5, tab6, tab7,
  fig3, fig4, fig5, fig6, fig7, fig8,
  file = here::here("experiments/tables_figures.RData")
)
