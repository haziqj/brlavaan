# Script to test SEM RBM functions using the HolzingerSwineford1939 dataset

library(tidyverse)
library(lavaan)
source("R/sem_rbm_functions.R")

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

# Check gradient
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

# Check Hessian
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

# Check cross product of scores
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
