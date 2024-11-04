library(tidyverse)
library(lavaan)
source(here::here("R", "10-sem_rbm_functions.R"))

HS.model <- "
  visual  =~ x1 + a*x2 + x3
  textual =~ x4 + a*x5 + x6
  speed   =~ x7 + a*x8 + x9
"

fit_lav   <- cfa(HS.model, HolzingerSwineford1939, information = "observed")
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

test_that("multiplication works", {
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
})
