skip()  # see test-growth2.R
## IK + OK, 29 Oct 2024
set.seed(111224)
dat <- gen_data_growth(n = 15, rel = 0.8, dist = "Normal", scale = 1 / 10)
mod <- txt_mod_growth(0.8)
tru <- truth(dat)

fit_lav    <- growth(mod, dat, start = tru)
fit_ML     <- fit_sem(mod, dat, start = tru, rbm = "none")
fit_eRBM   <- fit_sem(mod, dat, start = tru, rbm = "explicit")
fit_iRBM   <- fit_sem(mod, dat, start = tru, rbm = "implicit")
# fit_iRBM$converged

# fit_iRBMp  <- fit_sem(mod, dat, plugin_pen = pen_ridge)
# fit_iRBMpb <- fit_sem(mod, dat, plugin_pen = pen_ridge_bound)
# fit_huber  <- fit_sem(mod, dat, plugin_pen = pen_huber)

test_that("ML estimator matches lavaan", {
  expect_equal(
    as.numeric(coef(fit_lav)),
    coef(fit_ML),
    ignore_attr = TRUE,
    tolerance = 1e-4
  )
})

# tibble::tibble(
#   param = names(coef(fit_lav)),
#   truth = truth(dat),
#   ML = coef(fit_ML),
#   eRBM = coef(fit_eRBM),
#   iRBM = coef(fit_iRBM),
#   iRBMp_ridge = coef(fit_iRBMp),
#   iRBMp_ridgeb = coef(fit_iRBMpb),
#   iRBMp_huber = coef(fit_huber)
# )

N <- nrow(dat)
p <- ncol(dat)
Lambda <- matrix(c(rep(1, p), 1:p - 1), nrow = p, ncol = 2)

Loglik <- function(pars, dat, L, i = NULL) {
  q <- ncol(L)
  p <- nrow(L)
  N <- nrow(dat)
  alpha <- pars[1:q]
  theta <- pars[q + 1]
  psi <- pars[-(1:(q + 1))]
  sigma <- L %*% matrix(c(psi[1], psi[2], psi[2], psi[3]), q, q) %*% t(L) + diag(rep(theta, 10))
  ll_contr <- mvtnorm::dmvnorm(dat, mean = Lambda %*% alpha, sigma = sigma, log = TRUE)
  if (is.null(i)) {
    sum(ll_contr)
  } else {
    ll_contr[i]
  }
}

grad_Loglik <- function(pars, dat, L, i = NULL) {
  if (is.null(i)) {
    numDeriv::grad(Loglik, x = lav_c, dat = dat, L = Lambda, i = NULL)
  } else {
    numDeriv::grad(Loglik, x = lav_c, dat = dat, L = Lambda, i = i)
  }
}

emat <- function(pars, dat, L) {
  N <- nrow(dat)
  gr_contr <- sapply(1:N, function(i) grad_Loglik(pars, dat, L, i))
  tcrossprod(gr_contr)
}

jmat <- function(pars, dat, L) {
  -numDeriv::hessian(Loglik, x = lav_c, dat = dat, L = Lambda, i = NULL)
}

pen_Loglik <- function(pars, dat, L) {
  ll <- Loglik(pars, dat, L, i = NULL)
  e <- emat(pars, dat, L)
  jinv <- jmat(pars, dat, L) |> solve()
  ll - 0.5 * sum(diag(jinv %*% e))
}

# Getting the parameters in the correct order for the manual likelihood functions
ord_names <- c("i~1", "s~1", "v", "i~~i", "i~~s", "s~~s")
lav_c <- coef(fit_lav)[ord_names]

## ----- Start tests -----------------------------------------------------------

test_that("Check log-likelihood function against lavaan", {
  expect_equal(
    Loglik(lav_c, dat, Lambda),
    c(logLik(fit_lav))[1],
    ignore_attr = TRUE
  )
})

test_that("Check gradients against lavaan", {
  expect_equal(
    numDeriv::grad(Loglik, x = lav_c, dat = dat, L = Lambda),
    fit_lav@optim$dx,
    tolerance = 1e-05
  )
  expect_equal(
    numDeriv::grad(Loglik, x = lav_c, dat = dat, L = Lambda),
    grad_Loglik(lav_c, dat = dat, L = Lambda),
    tolerance = 1e-05
  )
})

test_that("Check Hessian against lavaan", {
  expect_equal(
    -numDeriv::hessian(Loglik, x = lav_c, dat = dat, L = Lambda),
    lavInspect(fit_lav, "hessian")[ord_names, ord_names] * lavaan::nobs(fit_lav),
    ignore_attr = TRUE,
    tolerance = 1e-01
  )
  expect_equal(
    -numDeriv::hessian(Loglik, x = lav_c, dat = dat, L = Lambda),
    jmat(lav_c, dat = dat, L = Lambda),
    ignore_attr = TRUE,
    tolerance = 1e-01
  )
})

# 15/2/2025: The following will not work... because brlavaan now uses
# Matrix::nearestPD() to fix the Hinv.

# test_that("Checking penalty", {
#   expect_equal(
#     brlavaan:::penalty(lav_c, fit_lav@Model, fit_lav@SampleStats, fit_lav@Data, fit_lav@Options, kind = "observed"),
#     pen_Loglik(lav_c, dat, Lambda) - Loglik(lav_c, dat, Lambda),
#     tolerance = 1e-4
#   )
# })

# neg_pen_Loglik <- function(pars, dat, L) -pen_Loglik(pars, dat, L)
# res1 <- nlminb(coef(fit_lav), neg_pen_Loglik, dat = dat, L = Lambda)
#
# test_that("Checking iRBM fit", {
#   expect_equal(res1$par, as.numeric(res2$par), ignore_attr = TRUE)
#   expect_equal(
#     res1$par,
#     coef(fit_iRBM),
#     ignore_attr = TRUE,
#     tolerance = 1e-02
#   )
# })

# test_that("Check iRBM fit", {
#   coef_pos <- c(2, 4, 6, 1, 5, 3)
#   expect_equal(
#     numDeriv::grad(pen_Loglik, coef(fit_iRBM)[coef_pos], dat = dat, L = Lambda),
#     rep(0, 6),
#     tolerance = 1e-02)
# })
