# Script to test old vs new functions

library(tidyverse)
library(lavaan)
library(tinytest)
source("R/10-sem_rbm_functions.R")
source("R/old/41-implicit_wrap.R")

n <- 1000
dat <- simulateData(
  model = "
    eta1 =~ 1*y1 + 0.8*y2 + 0.6*y3
    eta2 =~ 1*y4 + 0.8*y5 + 0.6*y6
    eta2 ~ 0.3*eta1
  ",
  sample.nobs = n
)

# {lavaan} fit
mod <- "
  eta1 =~ y1 + y2 + y3
  eta2 =~ y4 + y5 + y6
"
fit_lav   <- sem(mod, dat)
fit_ML    <- fit_sem(mod, dat, method = "ML")
fit_eRBM  <- fit_sem(mod, dat, method = "eRBM")
fit_iRBM  <- fit_sem(mod, dat, method = "iRBM")
fit_iRBMp <- fit_sem(mod, dat, method = "iRBMp")

m <- length(coef(fit_lav))
res_old <- imp_twofac(dat, trace = 10)

# Check old method agrees with {lavaan}
expect_equal(coef(fit_lav), res_old$ml, tol = 1e-4)
est_iRBMp_old <- res_old$iRBMp
class(est_iRBMp_old) <- "numeric"
expect_equal(coef(fit_iRBMp), est_iRBMp_old, tol = 1e-3)

# Check new method agrees with old method
est_ML_old <- res_old$ml
class(est_ML_old) <- "numeric"
expect_equal(coef(fit_ML), est_ML_old, tol = 1e-4)

