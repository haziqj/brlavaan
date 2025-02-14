skip()
# The iRBM method requires good starting values, especially for small sample
# sizes.

set.seed(456)
dat <- gen_data_twofac(n = 15, rel = 0.8, dist = "Normal")
mod <- txt_mod_twofac(0.8)

fit_ML     <- fit_sem(mod, dat, rbm = "none")
fit_eRBM   <- fit_sem(mod, dat, rbm = "explicit")
fit_iRBM   <- fit_sem(mod, dat, rbm = "implicit", info_se = "expected")

test_that("False convergence in iRBM", {expect_false(fit_iRBM$converged)})

test_that(
  "Convergence OK with good starting values",
  {
    fit_iRBM <- fit_sem(mod, dat, rbm = "implicit", start = coef(fit_ML))
    expect_true(fit_iRBM$converged)
    fit_iRBM <- fit_sem(mod, dat, rbm = "implicit", start = coef(fit_eRBM))
    expect_true(fit_iRBM$converged)
  }
)

# For the growth model, false convergences no matter what!
set.seed(789)
dat <- gen_data_growth(n = 15, rel = 0.8, dist = "Normal")
mod <- txt_mod_growth(0.8)

fit_ML     <- fit_sem(mod, dat, rbm = "none")
fit_eRBM   <- fit_sem(mod, dat, rbm = "explicit")
fit_iRBM   <- fit_sem(mod, dat, rbm = "implicit", info_se = "expected")

test_that("False convergence in iRBM", {expect_false(fit_iRBM$converged)})

test_that(
  "Convergence OK with good starting values",
  {
    fit_iRBM <- fit_sem(mod, dat, rbm = "implicit", start = coef(fit_ML))
    expect_false(fit_iRBM$converged)
    # expect_equal(coef(fit_iRBM), coef(fit_ML))

    fit_iRBM <- fit_sem(mod, dat, rbm = "implicit", start = coef(fit_eRBM))
    expect_false(fit_iRBM$converged)
    # expect_equal(coef(fit_iRBM), coef(fit_eRBM))
  }
)
