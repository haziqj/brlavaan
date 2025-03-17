dat <- gen_data_growth(n = 15, rel = 0.8, dist = "Normal", scale = 1 / 10)
mod <- txt_mod_growth(0.8)
tru <- truth(dat)

fit_lav   <- lavaan::growth(mod, dat, ceq.simple = TRUE)
fit_ML    <- fit_sem(mod, dat, rbm = "none", ceq.simple = TRUE)

list2env(get_lav_stuff(fit_lav), environment())

# Packed and unpacked version of coefficients
z.unpack <- coef(fit_lav)
z.pack   <- pack_x(z.unpack, lavmodel)

# Check log-likelihood
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
    loglik_val2 <- LOGLIK(z.pack, mod, dat)
    expect_equal(loglik_val, as.numeric(logLik(fit_lav)), tolerance = 1e-4)
    expect_equal(loglik_val, loglik_val2, tolerance = 1e-4)

  }
)

# Check scores
test_that(
  "Gradient is zero at the true values",
  {
    grad_lav <- grad_loglik(
      theta = z.pack,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
    grad_lav2 <- GRAD(z.pack, mod, dat)

    grad_num <- numDeriv::grad(
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
    grad_num2 <- numDeriv::grad(LOGLIK, z.pack, model = mod, data = dat)

    expect_equal(grad_lav, rep(0, length(grad_lav)), tolerance = 1e-4)
    expect_equal(grad_lav, grad_num, tolerance = 1e-4)
    expect_equal(grad_lav, grad_lav2, tolerance = 1e-4)
    expect_equal(grad_lav, grad_num2, tolerance = 1e-4)
    expect_equal(sign(grad_lav), sign(grad_num))
  }
)

test_that("Gradient function agrees", {
  x <- pack_x(tru, lavmodel)
  grad_fitsem <- grad_loglik(  # uses lav_model_gradient()
    theta = x,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
  grad_manual <- GRAD(x, mod, dat)  # uses closed form expressions

  numgrad_fitsem <- numDeriv::grad(
    func = loglik,
    x = x,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    bias_reduction = FALSE,
    kind = "observed",
    plugin_pen = NULL,
    verbose = FALSE
  )
  numgrad_manual <- numDeriv::grad(LOGLIK, x, model = mod, data = dat)

  expect_equal(grad_fitsem, grad_manual, tolerance = 1e-5)
  expect_equal(grad_fitsem, numgrad_fitsem, tolerance = 1e-5)
  expect_equal(grad_fitsem, numgrad_manual, tolerance = 1e-5)

  expect_equal(grad_manual, numgrad_manual, tolerance = 1e-5)
  expect_equal(grad_manual, numgrad_fitsem, tolerance = 1e-5)

  expect_equal(numgrad_fitsem, numgrad_manual, tolerance = 1e-5)

})
