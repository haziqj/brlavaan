test_that("brsem works", {
  mod <- txt_mod_twofac(0.8)
  dat <- gen_data_twofac(n = 50, rel = 0.8)
  fit <- brsem(mod, dat, estimator.args = list(rbm = FALSE))
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
               estimator.args = list(rbm = "iRBM", plugin_penalty = "pen_ridge"))
  tmp <- capture.output(print(fit))
  tmp <- capture.output(summary(fit))
  expect_true(is(fit, "brlavaan"))
  expect_true(is(coef(fit), "lavaan.vector"))
})

test_that("brgrowth works", {
  mod <- txt_mod_growth(0.5)
  dat <- gen_data_growth(n = 50)
  fit <- brgrowth(mod, dat, estimator.args = list(rbm = "eRBM"))
  tmp <- capture.output(print(fit))
  tmp <- capture.output(summary(fit))
  expect_true(is(fit, "brlavaan"))
  expect_true(is(coef(fit), "lavaan.vector"))
})

