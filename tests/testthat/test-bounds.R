n <- 15
rel <- 0.8
dist <- "Normal"
seed <- 123

dat <- gen_data_growth(n = n, rel = rel, dist = dist, seed = seed)
mod <- txt_mod_growth(rel = rel)
tru <- truth(dat)

test_that("No bounds", {
  fit_lav <- growth(mod, dat, bounds = "none")
  fit_ML  <- fit_sem(mod, dat, rbm = "none", bounds = "none")
  expect_equal(coef(fit_lav), fit_ML$coef, tolerance = 1e-5)
})

test_that("No bounds with starting values", {
  fit_lav <- growth(mod, dat, start = tru)
  fit_ML  <- fit_sem(mod, dat, rbm = "none", start = tru)
  expect_equal(coef(fit_lav), fit_ML$coef, tolerance = 1e-3)
})

test_that("Standard bounds", {
  fit_lav <- growth(mod, dat, bounds = "standard")
  fit_ML  <- fit_sem(mod, dat, rbm = "none", bounds = "standard")
  expect_equal(coef(fit_lav), fit_ML$coef, tolerance = 1e-3)
})

test_that("Standard bounds with starting values", {
  fit_lav <- growth(mod, dat, start = tru, bounds = "standard")
  fit_ML  <- fit_sem(mod, dat, rbm = "none", start = tru, bounds = "standard")
  expect_equal(coef(fit_lav), fit_ML$coef, tolerance = 1e-3)
})

test_that("ceq.simple = FALSE gets ignored", {
  fit_ML  <- fit_sem(mod, dat, rbm = "none", start = tru, bounds = "standard",
                     ceq.simple = FALSE)
  expect_true(length(fit_ML$bounds$lower) > 0)
  expect_true(length(fit_ML$bounds$upper) > 0)
})
