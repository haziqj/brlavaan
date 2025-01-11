# Script to test the main fit_sem() function using the HolzingerSwineford1939
# dataset. This mainly tests the maximum likelihood estimation.
library(lavaan)

HS.model <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
"
data(HolzingerSwineford1939, package = "lavaan")

fit_lav   <- lavaan::cfa(HS.model, HolzingerSwineford1939, information = "observed")
fit_ML    <- fit_sem(HS.model, HolzingerSwineford1939, rbm = "none")

list2env(get_lav_stuff(fit_lav), environment())

# Packed version of coefficients
z.unpack <- coef(fit_lav)
if (lavmodel@eq.constraints) {
  z.pack <- as.numeric(
    (z.unpack - lavmodel@eq.constraints.k0) %*% lavmodel@eq.constraints.K
  )
} else {
  z.pack <- z.unpack
}

# DEBUG ------------------------------------------------------------------------
# lavInspect(fit_lav, "free")
# tmp <- fit_sem(HS.model, HolzingerSwineford1939, debug = TRUE)

test_that(
  "Log-likelihood value matches lavaan::logLik()",
  {
    loglik_val <- loglik(
      theta = z.pack,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      bias_reduction = FALSE,
      kind = "observed",
      plugin_pen = NULL,
      verbose = FALSE
    )
    expect_equal(loglik_val, as.numeric(logLik(fit_lav)), tolerance = 1e-4)
  }
)

# Test maximum likelihood estimation -------------------------------------------
test_that(
  "Maximum likelihood estimation matches lavaan::sem()",
  {
    est_lav <- coef(fit_lav)
    class(est_lav) <- "numeric"
    est_ML <- coef(fit_ML)
    expect_equal(est_lav, est_ML, tolerance = 1e-4, ignore_attr = FALSE)

    sd_lav <- unname(sqrt(diag(vcov(fit_lav))))
    sd_ML <- fit_ML$stderr
    expect_equal(sd_lav, sd_ML, tolerance = 1e-4, ignore_attr = FALSE)
  }
)

# Test gradient ----------------------------------------------------------------
test_that(
  "Gradient matches numDeriv::grad()",
  {
    grad.lav <- grad_loglik(
      theta = z.pack,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
    grad.num <- numDeriv::grad(
      func = loglik,
      x = z.pack,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      bias_reduction = FALSE,
      kind = "observed",
      plugin_pen = NULL,
      verbose = FALSE
    )
    expect_equal(grad.lav, rep(0, length(grad.lav)), tolerance = 1e-4)
    expect_equal(grad.lav, grad.num, tolerance = 1e-4)
    expect_equal(sign(grad.lav), sign(grad.num))
  }
)

# Test Hessian -----------------------------------------------------------------
test_that(
  "Hessian matches numDeriv::hessian()",
  {
    hessian.lav <- hessian_loglik(
      theta = z.pack,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
    if (lavmodel@eq.constraints) {
      hessian.lav <- t(lavmodel@eq.constraints.K) %*% hessian.lav %*%
        lavmodel@eq.constraints.K
    }
    hessian.num1 <- numDeriv::hessian(
      func = loglik,
      x = z.pack,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      bias_reduction = FALSE,
      kind = "observed",
      plugin_pen = NULL,
      verbose = FALSE
    )
    hessian.num2 <- numDeriv::jacobian(
      func = grad_loglik,
      x = z.pack,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
    expect_equal(hessian.lav, hessian.num1, tolerance = 1e-5)
    expect_equal(hessian.lav, hessian.num2, tolerance = 1e-5)
  }
)

# Test cross product of scores -------------------------------------------------
test_that(
  "Cross product of scores matches lavaan:::crossprod(scores_loglik())",
  {
    scores.lav <- scores_loglik(
      theta = z.pack,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
    crossprod.lav <- crossprod(scores.lav)
    crossprod.num <- crossprod(scores_loglik(
      theta = z.pack,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    ))
    expect_equal(crossprod.lav, crossprod.num, tolerance = 1e-5)
  }
)
