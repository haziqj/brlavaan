library(brlavaan)
library(tidyverse)
library(furrr)
theme_set(theme_bw())
here::i_am("experiments/10-run_sims.R")
load(here::here("experiments/simu_id.RData"))

ncores <- future::availableCores() - 1
future::plan(multisession, workers = ncores)
B <- 1000  # Number of simulations

## ----- Run sims --------------------------------------------------------------
# source("experiments/11-run_sims_twofac.R")
source("experiments/12-run_sims_growth.R")
