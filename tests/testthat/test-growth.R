## IK + OK, 29 Oct 2024
set.seed(291024)
dat <- gen_data_growth(n = 100, rel = 0.5, dist = "Normal")
mod <- txt_mod_growth(0.5)

fit_lav   <- growth(mod, dat)
fit_ML    <- fit_sem(mod, dat, estimator = "ML", lavfun = "growth")
fit_eBRM  <- fit_sem(mod, dat, estimator = "eBRM", lavfun = "growth")
fit_iBRM  <- fit_sem(mod, dat, estimator = "iBRM", lavfun = "growth")
fit_iBRMp <- fit_sem(mod, dat, estimator = "iBRMp", lavfun = "growth")

test_that("ML estimator matches lavaan", {
  expect_equal(
    as.numeric(coef(fit_lav)),
    coef(fit_ML),
    ignore_attr = TRUE,
    tolerance = 1e-5
  )
})

# tibble(
#   param = names(coef(fit_lav)),
#   truth = truth(dat),
#   ML = coef(fit_ML),
#   eBRM = coef(fit_eBRM),
#   iBRM = coef(fit_iBRM),
#   iBRMp = coef(fit_iBRMp),
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
    tolerance = 1e-06
  )
  expect_equal(
    numDeriv::grad(Loglik, x = lav_c, dat = dat, L = Lambda),
    grad_Loglik(lav_c, dat = dat, L = Lambda),
    tolerance = 1e-06
  )
})

test_that("Check Hessian against lavaan", {
  expect_equal(
    -numDeriv::hessian(Loglik, x = lav_c, dat = dat, L = Lambda),
    lavInspect(fit_lav, "hessian")[ord_names, ord_names] * 100,
    ignore_attr = TRUE,
    tolerance = 1e-02
  )
  expect_equal(
    -numDeriv::hessian(Loglik, x = lav_c, dat = dat, L = Lambda),
    jmat(lav_c, dat = dat, L = Lambda),
    ignore_attr = TRUE,
    tolerance = 1e-02
  )
})

# neg_pen_Loglik <- function(pars, dat, L) -pen_Loglik(pars, dat, L)
# res1 <- nlminb(coef(fit_lav), neg_pen_Loglik, dat = dat, L = Lambda)
#
# test_that("Checking iBRM fit", {
#   expect_equal(res1$par, as.numeric(res2$par), ignore_attr = TRUE)
#   expect_equal(
#     res1$par,
#     coef(fit_iBRM),
#     ignore_attr = TRUE,
#     tolerance = 1e-02
#   )
# })

test_that("Check iBRM fit", {
  coef_pos <- c(2, 4, 6, 1, 5, 3)
  expect_equal(
    numDeriv::grad(pen_Loglik, coef(fit_iBRM)[coef_pos], dat = dat, L = Lambda),
    rep(0, 6),
    tolerance = 1e-02)
})
