library(tidyverse)
theme_set(theme_bw())
library(gt)
library(latex2exp)
here::i_am("experiments/20-analyse_sims.R")

mycols <- c(
  ML = "#E31A1C",
  eRBM = "#A6CEE3",
  iRBM = "#1F78B4",
  Jackknife = "#B2DF8A",
  Bootstrap = "#33A02C",
  `Ozenne et al.` = "#FDBF6F",
  REML = "#FF7F00"
)

## ----- Prepare results -------------------------------------------------------
source("experiments/21-prep_results.R")
# Creates:
# 1. res_nested [for easily finding out where convergence issues are]
# 2. res_full [raw results, combined]
# 3. res [raw results, not combined]
# 4. res_summ [all results, summarised (for plots)]
# 5. res_dr [results from Dhaene & Rosseel]
# 6. res [raw results from our simulations]
# 7. and others...

## ----- Analyse results -------------------------------------------------------
source(here::here("experiments/22-convergence.R"))
source(here::here("experiments/23-analysis_growth.R"))
source(here::here("experiments/24-analysis_twofac.R"))

## ----- Save ------------------------------------------------------------------
save(
  tab1, tab2, tab3, tab4, tab5, tab6, tab7,
  fig3, fig4, fig5, fig6, fig7, fig8,
  timing,
  fig_growth_dis, fig_growth_perf,
  fig_twofac_dis, fig_twofac_perf,
  file = here::here("experiments/tables_figures.RData")
)
