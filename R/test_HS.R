# Script to test SEM RBM functions using the HolzingerSwineford1939 dataset

library(tidyverse)
library(lavaan)
library(tinytest)
source("R/10-sem_rbm_functions.R")

HS.model <- "
  visual  =~ x1 + a*x2 + x3
  textual =~ x4 + a*x5 + x6
  speed   =~ x7 + a*x8 + x9
"

fit_lav   <- cfa(HS.model, HolzingerSwineford1939)
fit_ML    <- fit_sem(HS.model, HolzingerSwineford1939, method = "ML")
fit_eRBM  <- fit_sem(HS.model, HolzingerSwineford1939, method = "eRBM")
fit_iRBM  <- fit_sem(HS.model, HolzingerSwineford1939, method = "iRBM")
fit_iRBMp <- fit_sem(HS.model, HolzingerSwineford1939, method = "iRBMp")

# Packed version of coefficients
with(get_lav_stuff(fit_lav), {
  z.unpack <- coef(fit_lav)
  if (lavmodel@eq.constraints) {
    z.pack <<- as.numeric(
      (z.unpack - lavmodel@eq.constraints.k0) %*% lavmodel@eq.constraints.K
    )
  } else {
    z.pack <<- z.unpack
  }
})

# DEBUG ------------------------------------------------------------------------
lavInspect(fit_lav, "free")
tmp <- fit_sem(HS.model, HolzingerSwineford1939, debug = TRUE)

# Test likelihood value --------------------------------------------------------
with(get_lav_stuff(fit_lav), {
  loglik.val <<- loglik(
    theta = z.pack,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
})
expect_equal(loglik.val, as.numeric(logLik(fit_lav)), tol = 1e-4)

# Test maximum likelihood estimation -------------------------------------------
est_lav <- as.numeric(coef(fit_lav))
est_ML <- coef(fit_ML)
expect_equal(est_lav, est_ML, tol = 1e-4, check.attributes = FALSE)

sd_lav <- sqrt(diag(vcov(fit_lav)))
sd_ML <- fit_ML$stderr
expect_equal(sd_lav, sd_ML, tol = 1e-4, check.attributes = FALSE)

# Test gradient ----------------------------------------------------------------
with(get_lav_stuff(fit_lav), {
  grad.lav <<- grad_loglik(
    theta = z.pack,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
  grad.num <<- numDeriv::grad(
    func = loglik,
    x = z.pack,
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
    theta = z.pack,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
  hessian.num1 <<- -numDeriv::hessian(
    func = loglik,
    x = z.pack,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
  hessian.num2 <<- -numDeriv::jacobian(
    func = grad_loglik,
    x = z.pack,
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
    theta = z.pack,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
  e_mat2 <<- crossprod(scores_loglik(
    theta = z.pack,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  ))
})
expect_equal(e_mat1, e_mat2, tol = 1e-5)

# Take a look at all methods ---------------------------------------------------
tibble(
  id = rep(seq_along(coef(fit_lav)), 4),
  param = rep(names(coef(fit_lav)), 4),
  est = c(fit_ML$coef, fit_eRBM$coef, fit_iRBM$coef, fit_iRBMp$coef),
  se = c(fit_ML$stderr, fit_eRBM$stderr, fit_iRBM$stderr, fit_iRBMp$stderr),
  method = rep(c("ML", "eRBM", "iRBM", "iRBMp"), each = length(coef(fit_lav))),
  time = rep(c(fit_ML$time, fit_eRBM$time, fit_iRBM$time, fit_iRBMp$time),
             each = length(coef(fit_lav)))
) |>
  pivot_wider(id_cols = c(id, param), names_from = method,
              values_from = c(est, se)) |>
  select(-id) |>
  print(n = 100)

