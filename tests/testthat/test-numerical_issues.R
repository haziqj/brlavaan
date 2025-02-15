test_that("Near PD may fix convergences", {
  set.seed(21)
  dat <- gen_data_twofac(n = 15, rel = 0.5, dist = "Normal", scale = 1)
  mod <- txt_mod_twofac(0.5)
  tru <- truth(dat)

  fit1 <- suppressWarnings(fit_sem(mod, dat, rbm = "implicit", start = tru, nearPD = FALSE))
  expect_false(fit1$converged)
  fit2 <- fit_sem(mod, dat, rbm = "implicit", start = tru, nearPD = TRUE)
  expect_true(fit2$converged)
})

test_that("Data scaling for growth model", {
  # For the growth model, Interesting to see that DATA SCALING fixes the issue!

  mod <- txt_mod_growth(0.8)
  dat <- gen_data_growth(n = 15, rel = 0.8, dist = "Normal", seed = 789)
  fit1 <- fit_sem(mod, dat, start = truth(dat), rbm = "implicit")
  dat <- gen_data_growth(n = 15, rel = 0.8, dist = "Normal", seed = 789,
                         scale = 1/10)
  fit2 <- fit_sem(mod, dat, start = truth(dat), rbm = "implicit")

  expect_false(fit1$converged)
  expect_true(fit2$converged)
})
