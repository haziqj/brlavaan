# Script to test Growth Curve Model functions

library(tidyverse)
library(lavaan)
library(tinytest)
source("R/10-sem_rbm_functions.R")
source("R/20-gen_data.R")

# Reliability == 0.8 -----------------------------------------------------------
dat <- gen_data_growth(n = 15, rel = 0.8, dist = "Normal")
mod <- txt_mod_growth(0.8)
fit_lav   <- sem(mod, dat)

# Test contraction and expansion of theta
theta_full <- coef(fit_lav)
class(theta_full) <- "numeric"
pt <- partable(fit_lav)
idx <- get_constr_idx(pt)
theta_short <- contract_theta(theta_full, idx)
theta_full2 <- expand_theta(theta_short, idx)
expect_equal(theta_short, theta_full[unique(names(theta_full))])
expect_equal(theta_full, theta_full2)

fit_ML    <- fit_sem(mod, dat, method = "ML")
fit_eRBM  <- fit_sem(mod, dat, method = "eRBM")
fit_iRBM  <- fit_sem(mod, dat, method = "iRBM")
fit_iRBMp <- fit_sem(mod, dat, method = "iRBMp")

# Reliability == 0.5 -----------------------------------------------------------
dat <- gen_data_growth(n = 1000, rel = 0.5, dist = "Non-normal")
mod <- txt_mod_growth(0.5)

fit_lav   <- sem(mod, dat)
fit_ML    <- fit_sem(mod, dat, method = "ML")
fit_eRBM  <- fit_sem(mod, dat, method = "eRBM")
fit_iRBM  <- fit_sem(mod, dat, method = "iRBM")
fit_iRBMp <- fit_sem(mod, dat, method = "iRBMp")

tibble(
  param = names(coef(fit_lav)),
  truth = truth(dat),
  ML = coef(fit_ML),
  eRBM = coef(fit_eRBM),
  iRBM = coef(fit_iRBM),
  iRBMp = coef(fit_iRBMp),
)
