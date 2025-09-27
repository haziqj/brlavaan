set.seed(1)
n <- 15
rel <- 0.8
dist <- "Normal"
dat <- gen_data_twofac(n = n, rel = rel, dist = dist, scale = 1 / 10)
mod <- txt_mod_twofac(rel)
tru <- truth(dat)

fit_lav <- lavaan::sem(mod, dat)
fit_ML  <- fit_sem(mod, dat, rbm = "none")
fit_eBR <- fit_sem(mod, dat, rbm = "explicit", start = tru)
fit_iBR <- fit_sem(mod, dat, rbm = "implicit", start = tru)

# list(fit_ML = fit_ML, fit_eBR = fit_eBR, fit_iBR = fit_iBR, fit_lav = fit_lav) |>
#   map_dbl(\(x) (coef(x)["i~~s"] - tru["i~~s"]) / tru["i~~s"])

# Check log-likelihood
test_that("Log-likelihood value matches lavaan::logLik()", {
  with(get_lav_stuff(fit_lav), {
    x <<- lavaan::lav_model_get_parameters(lavmodel)
    loglik_val <<- brlavaan:::loglik(
      x = x,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      bias_reduction = FALSE,
      plugin_pen = NULL,
      verbose = FALSE
    )
  })
  loglik_val2 <- brlavaan:::LOGLIK_TWOFAC(x, mod, dat)

  expect_equal(loglik_val, as.numeric(logLik(fit_lav)), tolerance = 1e-4)
  expect_equal(loglik_val, loglik_val2, tolerance = 1e-4)
})

# Check scores are zero at optima
test_that("Gradient is zero at the true values", {

  # brlavaan internals
  with (get_lav_stuff(fit_lav), {
    x <<- lav_model_get_parameters(lavmodel)
    grad_lav <<- brlavaan:::grad_loglik(  # uses lav_model_gradient()
      x = x,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
    grad_num <<- numDeriv::grad(
      func = brlavaan:::loglik,
      x = x[!is.na(x)],
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      bias_reduction = FALSE,
      plugin_pen = NULL,
      bounds = bounds,
      verbose = FALSE
    )
  })

  # manual calculations
  grad_lav2 <- brlavaan:::GRAD_TWOFAC(x, mod, dat)
  grad_num2 <- numDeriv::grad(brlavaan:::LOGLIK_TWOFAC, x, model = mod, data = dat)

  expect_equal(grad_lav, rep(0, length(grad_lav)), tolerance = 1e-2)
  expect_equal(grad_lav, grad_num, tolerance = 1e-3)

  expect_equal(grad_lav, grad_lav2, tolerance = 1e-4)
  expect_equal(grad_lav, grad_num2, tolerance = 1e-2)

  expect_equal(sign(grad_lav), sign(grad_num))
})

test_that("Gradient function agrees", {

  x_lavmodel <- lavaan::lav_model_set_parameters(fit_lav@Model, x = c(tru, rep(0, 6)))
  x <- lavaan::lav_model_get_parameters(x_lavmodel)

  # brlavaan internals
  with (get_lav_stuff(fit_lav), {
    grad_lav <<- brlavaan:::grad_loglik(  # uses lav_model_gradient()
      x = x,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
    grad_num <<- numDeriv::grad(
      func = brlavaan:::loglik,
      x = x,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      bias_reduction = FALSE,
      plugin_pen = NULL,
      bounds = bounds,
      verbose = FALSE
    )
  })

  # manual calculations
  grad_lav2 <- brlavaan:::GRAD_TWOFAC(x[1:13], mod, dat)
  grad_num2 <- numDeriv::grad(brlavaan:::LOGLIK_TWOFAC, x[1:13], model = mod, data = dat)

  expect_equal(grad_lav, grad_num, tolerance = 1e-5)
  expect_equal(grad_lav, grad_lav2, tolerance = 1e-5)
  expect_equal(grad_lav, grad_num2, tolerance = 1e-5)

  expect_equal(grad_lav2, grad_num2, tolerance = 1e-5)
  expect_equal(grad_lav2, grad_num, tolerance = 1e-5)

  expect_equal(grad_num, grad_num2, tolerance = 1e-5)

})

# Check E matrix
test_that("E matrices agree", {

  x_lavmodel <- lavaan::lav_model_set_parameters(fit_lav@Model,
                                                 x = c(tru, rep(0, 6)))
  x <- lavaan::lav_model_get_parameters(x_lavmodel)

  e1 <- with(get_lav_stuff(fit_lav), {
    brlavaan:::information_matrix(  # uses lav_model_information_firstorder()
      x = x,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      kind = "first.order"
    )
  })
  e2 <- brlavaan:::EMAT_TWOFAC(x, mod, dat)

  expect_equal(e1, e2, tolerance = 1e-5)
})

# Check J matrix
test_that("J matrices agree", {

  x_lavmodel <- lavaan::lav_model_set_parameters(fit_lav@Model, x = tru)
  x <- lavaan::lav_model_get_parameters(x_lavmodel)

  j1 <- with(get_lav_stuff(fit_lav), {
    brlavaan:::information_matrix(  # uses lav_model_information_observed()
      x = x,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      kind = "observed"
    )
  })

  j2 <- brlavaan:::JMAT_TWOFAC(x, mod, dat)

  expect_equal(j1, j2, tolerance = 1e-4)
})

# Check penalty term
test_that("Penalty term correct", {

  pen1 <- with(get_lav_stuff(fit_lav), {
    x <<- lavaan::lav_model_get_parameters(lavmodel)

    brlavaan:::penalty(
      x = x,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
  })

  pen2 <- brlavaan:::PENALTY_TWOFAC(x, mod, dat)

  expect_equal(pen1, pen2, tolerance = 1e-2)

})

# Check explicit bias correction
test_that("eBR bias correction is correct", {

  xx <- coef(fit_lav)
  x_lavmodel <- lavaan::lav_model_set_parameters(fit_lav@Model, x = xx)
  x <- lavaan::lav_model_get_parameters(x_lavmodel)

  bias1 <- with(get_lav_stuff(fit_lav), {
    brlavaan:::bias(
      x = x,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
  })
  bias2 <- brlavaan:::BIAS_TWOFAC(x, mod, dat)

  expect_equal(bias1, bias2, tolerance = 1e-2)
})
