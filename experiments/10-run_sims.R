library(brlavaan)
library(tidyverse)
library(furrr)
theme_set(theme_bw())

ncores <- future::availableCores() - 1
future::plan(multisession, workers = ncores)
B <- 10  # Number of simulations

## ----- Run sims --------------------------------------------------------------
source("R/11-run_sims_twofac.R")
source("R/12-run_sims_growth.R")
