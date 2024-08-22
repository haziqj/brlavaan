# Script to test Growth Curve Model functions

library(tidyverse)
library(lavaan)
library(tinytest)
source("R/10-sem_rbm_functions.R")
source("R/20-gen_data.R")

# Reliability == 0.8 -----------------------------------------------------------
dat <- gen_data_growth(n = 100, rel = 0.8, dist = "Normal")
mod <- txt_mod_growth(0.8)

fit_lav   <- sem(mod, dat)
fit_ML    <- fit_sem(mod, dat, method = "ML")
fit_eRBM  <- fit_sem(mod, dat, method = "eRBM")
fit_iRBM  <- fit_sem(mod, dat, method = "iRBM")
fit_iRBMp <- fit_sem(mod, dat, method = "iRBMp")
