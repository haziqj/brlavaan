library(brlavaan)
library(tidyverse)
library(furrr)
theme_set(theme_bw())
load(here::here("experiments/simu_id.RData"))

ncores <- future::availableCores() - 1
future::plan(multisession, workers = ncores)
B <- 2000  # Number of simulations
