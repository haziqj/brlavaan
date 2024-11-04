library(tidyverse)
library(lavaan)
library(furrr)
theme_set(theme_bw())
source("R/12-sem_rbm_functions.R")
source("R/20-gen_data.R")
source("R/21-sim_functions.R")

ncores <- future::detectCores() - 2
future::plan(multisession, workers = ncores)
nsimu <- 1000

## ----- Run sims --------------------------------------------------------------
source("R/31-run_sims_twofac.R")
source("R/32-run_sims_growth.R")
