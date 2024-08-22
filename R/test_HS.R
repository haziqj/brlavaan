# Script to test SEM RBM functions using the HolzingerSwineford1939 dataset

library(tidyverse)
library(lavaan)
library(tinytest)
source("R/10-sem_rbm_functions.R")

HS.model <- "
  visual  =~ x1 + x2 + x3
  # textual =~ x4 + x5 + x6
  # speed   =~ x7 + x8 + x9
"

fit_lav   <- cfa(HS.model, HolzingerSwineford1939)
fit_ML    <- fit_sem(HS.model, HolzingerSwineford1939, method = "ML")
fit_eRBM  <- fit_sem(HS.model, HolzingerSwineford1939, method = "eRBM")
fit_iRBM  <- fit_sem(HS.model, HolzingerSwineford1939, method = "iRBM")
fit_iRBMp <- fit_sem(HS.model, HolzingerSwineford1939, method = "iRBMp")

# Test maximum likelihood estimation -------------------------------------------
est_lav <- coef(fit_lav)
class(est_lav) <- "numeric"
est_ML <- coef(fit_ML)
expect_equal(est_lav, est_ML, tol = 1e-4)

sd_lav <- sqrt(diag(vcov(fit_lav)))
sd_ML <- fit_ML$stderr
expect_equal(sd_lav, sd_ML, tol = 1e-4, check.attributes = FALSE)

# Test gradient ----------------------------------------------------------------
with(get_lav_stuff(fit_lav), {
  grad.lav <<- grad_loglik(
    theta = coef(fit_lav),
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
  grad.num <<- numDeriv::grad(
    func = loglik,
    x = coef(fit_lav),
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
})
expect_equal(grad.lav, rep(0, length(grad.lav)), tol = 1e-4)
expect_equal(grad.lav, grad.num, tol = 1e-4)
expect_equal(sign(grad.lav), sign(grad.num))

# Test Hessian -----------------------------------------------------------------
with(get_lav_stuff(fit_lav), {
  hessian.lav <<- hessian_loglik(
    theta = coef(fit_lav),
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
  hessian.num1 <<- -numDeriv::hessian(
    func = loglik,
    x = coef(fit_lav),
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
  hessian.num2 <<- -numDeriv::jacobian(
    func = grad_loglik,
    x = coef(fit_lav),
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
})
expect_equal(hessian.lav, hessian.num1, tol = 1e-5)
expect_equal(hessian.lav, hessian.num2, tol = 1e-5)

# Test cross product of scores -------------------------------------------------
with(get_lav_stuff(fit_lav), {
  e_mat1 <<- first_order_unit_information_loglik(
    theta = coef(fit_lav),
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
  e_mat2 <<- crossprod(scores_loglik(
    theta = coef(fit_lav),
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  ))
})
expect_equal(e_mat1, e_mat2, tol = 1e-5)
