test_that("Data scaling for growth model", {
  # For the growth model, Interesting to see that DATA SCALING fixes the issue!
  # 18/3/24 update: Able to use bounds = "standard" to stabilise optimisation

  mod <- txt_mod_growth(0.8)
  dat <- gen_data_growth(n = 15, rel = 0.8, dist = "Normal", seed = 789)
  fit1 <- fit_sem(mod, dat, start = truth(dat), rbm = "implicit")
  fit1b <- fit_sem(mod, dat, start = truth(dat), rbm = "implicit", bounds = "standard")
  dat <- gen_data_growth(n = 15, rel = 0.8, dist = "Normal", seed = 789,
                         scale = 1/10)
  fit2 <- fit_sem(mod, dat, start = truth(dat), rbm = "implicit")

  expect_false(fit1$converged)
  expect_true(fit1b$converged)
  expect_true(fit2$converged)
})
