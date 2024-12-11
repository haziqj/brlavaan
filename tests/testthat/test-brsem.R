test_that("brsem works", {
  mod <- txt_mod_twofac(0.8)
  dat <- gen_data_twofac(n = 50, rel = 0.8)
  fit <- brsem(mod, dat, estimator.args = list(rbm = "none"))
  expect_true(is(fit, "brlavaan"))
  expect_true(is(coef(fit), "lavaan.vector"))
  expect_snapshot(summary(fit))
})

test_that("brcfa works", {
  HS.model <- "
    visual  =~ x1 + x2 + x3
    # textual =~ x4 + x5 + x6
    # speed   =~ x7 + x8 + x9
  "
  fit <- brcfa(HS.model, data = HolzingerSwineford1939,
               estimator.args = list(rbm = "implicit", plugin_penalty = pen_ridge))

  expect_no_error({tmp <- capture.output(print(fit))})
  expect_no_error({tmp <- capture.output(summary(fit))})
  expect_true(is(fit, "brlavaan"))
  expect_true(is(coef(fit), "lavaan.vector"))
})

test_that("brgrowth works", {
  mod <- txt_mod_growth(0.5)
  dat <- gen_data_growth(n = 50)
  fit <- brgrowth(mod, dat, estimator.args = list(rbm = "explicit"))

  expect_no_error({tmp <- capture.output(print(fit))})
  expect_no_error({tmp <- capture.output(summary(fit))})
  expect_true(is(fit, "brlavaan"))
  expect_true(is(coef(fit), "lavaan.vector"))
})

test_that("brsem works", {
  HS.model <- "
    visual  =~ x1 + x2 + x3
    # textual =~ x4 + x5 + x6
    # speed   =~ x7 + x8 + x9
  "
  ML <- list(rbm = FALSE)  # should be "none" but this is a test
  eRBM  <- list(rbm = "explicit")
  iRBM  <- list(rbm = "implicit")
  iRBMp <- list(rbm = "implicit", plugin_penalty = pen_ridge)

  fit_ML <- brcfa(HS.model, data = HolzingerSwineford1939, estimator.args = ML)
  fit_eRBM <- brcfa(HS.model, data = HolzingerSwineford1939, estimator.args = eRBM)
  fit_iRBM <- brcfa(HS.model, data = HolzingerSwineford1939, estimator.args = iRBM)
  fit_iRBMp <- brcfa(HS.model, data = HolzingerSwineford1939, estimator.args = iRBMp)

  expect_true(is_ML(fit_ML))
  expect_true(is_eRBM(fit_eRBM))
  expect_true(is_iRBM(fit_iRBM))
  expect_true(is_iRBMp(fit_iRBMp, quietly = TRUE))
})
