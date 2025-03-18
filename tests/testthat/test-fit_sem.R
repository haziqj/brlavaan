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

test_that("Log-likelihood value matches lavaan::logLik()", {
  loglik_val <- with(get_lav_stuff(fit_lav), {
    loglik(
      x = coef(fit_lav),
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      bias_reduction = FALSE,
      plugin_pen = NULL,
      verbose = FALSE
    )
  })
  expect_equal(loglik_val, as.numeric(logLik(fit_lav)), tolerance = 1e-4)
})

test_that("Maximum likelihood estimation matches lavaan::sem()", {
  est_lav <- coef(fit_lav)
  class(est_lav) <- "numeric"
  est_ML <- coef(fit_ML)
  expect_equal(est_lav, est_ML, tolerance = 1e-4, ignore_attr = FALSE)

  sd_lav <- unname(sqrt(diag(vcov(fit_lav))))
  sd_ML <- fit_ML$stderr
  expect_equal(sd_lav, sd_ML, tolerance = 1e-4, ignore_attr = FALSE)
})

test_that("Gradient matches numDeriv::grad()", {
  with(get_lav_stuff(fit_lav), {
    grad_lav <<- grad_loglik(
      x = coef(fit_lav),
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
    grad_num <<- numDeriv::grad(
      func = loglik,
      x = coef(fit_lav),
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      bias_reduction = FALSE,
      plugin_pen = NULL,
      verbose = FALSE
    )
  })

  expect_equal(grad_lav, rep(0, length(grad_lav)), tolerance = 1e-4)
  expect_equal(grad_lav, grad_num, tolerance = 1e-4)
  expect_equal(sign(grad_lav), sign(grad_num))
})

test_that("Hessian matches numDeriv::hessian()", {
  with(get_lav_stuff(fit_lav), {
    hessian_lav <<- -1 * information_matrix(
      x = coef(fit_lav),
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      kind = "observed"
    )
    hessian_num1 <<- numDeriv::hessian(
      func = loglik,
      x = coef(fit_lav),
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      bias_reduction = FALSE,
      plugin_pen = NULL,
      verbose = FALSE
    )
    hessian_num2 <<- numDeriv::jacobian(
      func = grad_loglik,
      x = coef(fit_lav),
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
  })

  expect_equal(hessian_lav, hessian_num1, tolerance = 1e-5)
  expect_equal(hessian_lav, hessian_num2, tolerance = 1e-5)
})

test_that("Cross product of scores", {
  with(get_lav_stuff(fit_lav), {
    x <- lavaan::lav_model_get_parameters(lavmodel)
    lavimplied <- lavaan::lav_model_implied(lavmodel)
    moments <- list(cov = lavimplied$cov[[1]])
    if(lavmodel@meanstructure) moments$mean <- lavimplied$mean[[1]]
    ntab <- unlist(lavdata@norig)
    ntot <- sum(ntab)
    npar <- length(x)
    E1 <- lav_scores_ml(
      ntab = ntab,
      ntot = ntot,
      npar = npar,
      moments = moments,
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      lavmodel = lavmodel,
      lavoptions = lavoptions,
      scaling = FALSE
    )
    E1 <<- crossprod(E1)
    E2 <<- information_matrix(
      x = x,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      kind = "first.order"
    )
  })

  expect_equal(E1, E2, tolerance = 1e-5)
})
