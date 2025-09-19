set.seed(17)
N <- 100
REL <- 0.8
MOD <- "twofac"

if (MOD == "twofac") {
  dat <- gen_data_twofac(n = N, rel = REL, dist = "Normal")
  mod <- txt_mod_twofac(REL)
} else if (MOD== "growth") {
  dat <- gen_data_growth(n = N, rel = REL, dist = "Normal")
  mod <- txt_mod_growth(REL)
}
tru <- truth(dat)

test_that("No standard errors", {
  fit <- fit_sem(mod, dat, rbm = "none", se = "none")
  expect_true(all(is.na(fit$stderr)))
  expect_null(fit$vcov)
})

test_that("Standard argument for se", {
  # Observed information
  fit1 <- sem(mod, dat, se = "standard", information = "observed")
  fit2 <- fit_sem(mod, dat, rbm = "none", se = "standard", information = "observed", debug = !TRUE)
  fit3 <- brsem(mod, dat, estimator.args = list(rbm = "none"),
                se = "standard")  # default information is observed
  expect_equal(sqrt(diag(fit1@vcov$vcov)), fit2$stderr, tolerance = 1e-5)
  expect_equal(sqrt(diag(fit3@vcov$vcov)), fit2$stderr, tolerance = 1e-5)

  # Expected information
  fit1 <- sem(mod, dat, se = "standard", information = "expected")
  fit2 <- fit_sem(mod, dat, rbm = "none", se = "standard", information = "expected")
  fit3 <- brsem(mod, dat, estimator.args = list(rbm = "none"), information = "expected",
                se = "standard")
  expect_equal(sqrt(diag(fit1@vcov$vcov)), fit2$stderr, tolerance = 1e-5)
  expect_equal(sqrt(diag(fit3@vcov$vcov)), fit2$stderr, tolerance = 1e-5)

})

test_that("Robust HW argument for se", {
  # Observed information
  fit1 <- sem(mod, dat, se = "robust.huber.white", information = "observed")
  fit2 <- fit_sem(mod, dat, rbm = "none", se = "robust.huber.white", information = "observed")
  fit3 <- brsem(mod, dat, estimator.args = list(rbm = "none"), information = "observed",
                se = "robust.huber.white")  # default information is observed
  expect_equal(sqrt(diag(fit1@vcov$vcov)), fit2$stderr, tolerance = 1e-5)
  expect_equal(sqrt(diag(fit3@vcov$vcov)), fit2$stderr, tolerance = 1e-5)

  # Expected information
  fit1 <- sem(mod, dat, se = "robust.huber.white", information = "expected")
  fit2 <- fit_sem(mod, dat, rbm = "none", se = "robust.huber.white", information = "expected")
  fit3 <- brsem(mod, dat, estimator.args = list(rbm = "none"), information = "expected",
                se = "robust.huber.white")  # default information is observed
  expect_equal(sqrt(diag(fit1@vcov$vcov)), fit2$stderr, tolerance = 1e-5)
  expect_equal(sqrt(diag(fit3@vcov$vcov)), fit2$stderr, tolerance = 1e-5)

})

test_that("Other arguments for se should fail", {
  expect_error(fit_sem(mod, dat, se = "robust.sem"))
  expect_error(fit_sem(mod, dat, se = "bootstrap"))
})

