library(tidyverse)
theme_set(theme_bw())
library(gt)
library(latex2exp)
here::i_am("experiments/20-analyse_sims.R")
load(here::here("experiments/simu_res_growth.RData"))
load(here::here("experiments/simu_res_twofac.RData"))
mycols <- c(
  ML = "#E31A1C",
  eRBM = "#A6CEE3",
  iRBM = "#1F78B4",
  Jackknife = "#B2DF8A",
  Bootstrap = "#33A02C",
  `Ozenne et al.` = "#FDBF6F",
  REML = "#FF7F00"
)

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
    method = factor(method, levels = c("ML", "eRBM", "iRBM")),
    dist = factor(dist, levels = c("Normal", "Kurtosis", "Non-normal")),
    rel = factor(rel, levels = c("0.8", "0.5"), labels = c("Rel = 0.8", "Rel = 0.5")),
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

## ----- Growth & Two-factor results -------------------------------------------
source(here::here("experiments/21-analysis_growth.R"))
source(here::here("experiments/22-analysis_twofac.R"))

## ----- Save ------------------------------------------------------------------
save(
  tab1, tab2, tab3, tab4, tab5, tab6, tab7,
  fig3, fig4, fig5, fig6, fig7, fig8,
  timing,
  fig_growth_dis, fig_growth_perf,
  fig_twofac_dis, fig_twofac_perf,
  file = here::here("experiments/tables_figures.RData")
)
