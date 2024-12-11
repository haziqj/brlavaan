library(brlavaan)
library(tidyverse)
library(furrr)
theme_set(theme_bw())

ncores <- future::availableCores() - 1
future::plan(multisession, workers = ncores)
B <- 1000  # Number of simulations

## ----- Run sims --------------------------------------------------------------
source("experiments/11-run_sims_twofac.R")
source("experiments/12-run_sims_growth.R")
