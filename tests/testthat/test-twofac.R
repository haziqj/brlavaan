# Written by Ollie Kemp
set.seed(123)
dat <- gen_data_twofac(n = 15, rel = 0.8, dist = "Normal")
mod <- txt_mod_twofac(0.8)

eRBM  <- list(rbm = "eRBM")
iRBM  <- list(rbm = "iRBM")
iRBMp <- list(rbm = "iRBM", plugin_penalty = pen_ridge)

fit_lav   <- sem(mod, dat)
fit_ML    <- fit_sem(mod, dat)
fit_eRBM  <- fit_sem(mod, dat, estimator.args = eRBM, information = "observed")
fit_iRBM  <- fit_sem(mod, dat, estimator.args = iRBM, information = "expected")
fit_iRBMp <- fit_sem(mod, dat,  estimator.args = iRBMp, information = "expected")

N <- nrow(dat)
p <- ncol(dat)

# For this likelihood, we input the parameters 1:13 in the same order as lavaan
# output
Loglik <- function(pars, dat, i = NULL) {
  p <- ncol(dat)
  N <- nrow(dat)

  nu <- colMeans(dat)

  L <- matrix(c(1, pars[1], pars[2],0,0,0,0,0,0, 1, pars[3], pars[4]), nrow = p, ncol = 2)
  B <-  matrix(c(0, pars[5], 0, 0), nrow = 2, ncol = 2)
  theta <- diag(pars[6:11])
  psi <- diag(pars[12:13])
  IB_inv <- solve(diag(2) - B)

  sigma <- (L %*% IB_inv %*% psi %*% t(IB_inv) %*% t(L)) + theta

  ll_contr <- mvtnorm::dmvnorm(dat, mean=nu, sigma = sigma, log = TRUE)

  # Fully manual matches using mvtnorm
  #ll <- -0.5*N*log(det(sigma)) - N*p*0.5*log(2*pi)
  #for (j in 1:N) {
  #  ll <- ll - 0.5*as.matrix(dat[j, ]) %*% solve(sigma) %*% t(as.matrix(dat[j, ]))
  #}
  #print(ll)
  #print(unname(ll_contr))

  if (is.null(i)) {
    sum(ll_contr)
  } else {
    ll_contr[i]
  }
}

grad_Loglik <- function(pars, dat, i = NULL) {
  if (is.null(i)) {
    numDeriv::grad(Loglik, x = lav_c, dat = dat, i = NULL)
  } else {
    numDeriv::grad(Loglik, x = lav_c, dat = dat,  i = i)
  }
}

emat <- function(pars, dat) {
  N <- nrow(dat)
  gr_contr <- sapply(1:N, function(i) grad_Loglik(pars, dat, i))
  tcrossprod(gr_contr)
}

jmat <- function(pars, dat) {
  -numDeriv::hessian(Loglik, x = lav_c, dat = dat, i = NULL)
}

pen_Loglik <- function(pars, dat) {
  ll <- Loglik(pars, dat, i = NULL)
  e <- emat(pars, dat)
  jinv <- jmat(pars, dat) |> solve()
  ll - 0.5 * sum(diag(jinv %*% e))
}

neg_pen_Loglik <- function(pars, dat) {
  -pen_Loglik(pars, dat)
}

ord_names <- names(coef(fit_lav))
lav_c <- coef(fit_lav)[ord_names]
coef_pos <- c(1:13)

## ----- Start checks ----------------------------------------------------------
test_that("Checking likelihoods", {
  expect_equal(
    Loglik(lav_c, dat),
    c(logLik(fit_lav))[1],
    ignore_attr = FALSE
  )
  expect_equal(
    loglik(lav_c, fit_lav@Model, fit_lav@SampleStats, fit_lav@Data,
           fit_lav@Options, plugin_penalty = NULL, bias_reduction = FALSE,
           verbose = FALSE),
    Loglik(lav_c, dat)
  )
})

test_that("Checking gradients", {
  expect_equal(
    numDeriv::grad(Loglik, x = lav_c, dat = dat),
    fit_lav@optim$dx,
    tolerance = 1e-04
  )
  expect_equal(
    numDeriv::grad(Loglik, x = lav_c, dat = dat),
    grad_Loglik(lav_c, dat = dat),
    tolerance = 1e-04
  )
  expect_equal(
    grad_loglik(lav_c, fit_lav@Model, fit_lav@SampleStats, fit_lav@Data, fit_lav@Options),
    grad_Loglik(lav_c, dat),
    tolerance=1e-04
  )
})

test_that("Checking Hessian", {
  N <- nrow(dat)

  manual_hess <- -numDeriv::hessian(Loglik, x = lav_c, dat = dat)

  expect_equal(
    manual_hess,
    lavInspect(fit_lav, "hessian")[ord_names, ord_names] * N,
    ignore_attr = TRUE,
    tolerance = 1e-06
  )
  expect_equal(
    manual_hess,
    jmat(lav_c, dat = dat),
    ignore_attr = TRUE,
    tolerance = 1e-06
  )
  expect_equal(
    -manual_hess,
    hessian_loglik(lav_c, fit_lav@Model, fit_lav@SampleStats, fit_lav@Data, fit_lav@Options),
    ignore_attr = TRUE,
    tolerance = 1e-06
  )
})

test_that("Checking penalty", {
  expect_equal(
    penalty(lav_c, fit_lav@Model, fit_lav@SampleStats, fit_lav@Data, fit_lav@Options),
    pen_Loglik(lav_c, dat) - Loglik(lav_c, dat)
  )
})

# Checking the iRBM fit against the same using our implementation
res1 <- nlminb(coef(fit_lav), neg_pen_Loglik, dat = dat)
res2 <- optim(coef(fit_lav), pen_Loglik, dat = dat, method = "BFGS",
              control = list(fnscale = -1))

test_that("Checking iRBM fit", {
  expect_equal(res1$par, as.numeric(res2$par), ignore_attr = TRUE)
  expect_equal(
    res1$par,
    coef(fit_iRBM),
    ignore_attr = TRUE,
    tolerance = 1e-01
  )
})

# # Check gradients zero at optima
# expect_equal(
#   numDeriv::grad(pen_Loglik, coef(fit_iRBM)[coef_pos], dat = dat),
#   rep(0, 13),
#   tolerance = 1e-03
# )
# expect_equal(
#   numDeriv::grad(loglik, coef(fit_iRBM), lavmodel = fit_lav@Model, lavsamplestats = fit_lav@SampleStats, lavdata = fit_lav@Data, lavoptions = fit_lav@Options, RBM = TRUE),
#   rep(0, 13),
#   tolerance = 1e-03
# )
#
# # These parameters seem correct, at least wrt my penalised log likelihood
# expect_equal(numDeriv::grad(pen_Loglik, res1$par, dat = dat),
#              rep(0, 13),
#              tolerance = 1e-03)
