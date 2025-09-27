# Adapting the manual script for the growth model by HJ to the factor model
# Manually obtaining derivatives for the two factor model to check correctness

make_matrices_twofac <- function(x, uselavaan = FALSE, model, data) {
  # model and data arguments used only for creating fit0
  if (isTRUE(uselavaan)) {
    fit0 <- lavaan::sem(model = model, data = data, do.fit = FALSE,
                        meanstructure = FALSE)
    lavmodel <- lavaan::lav_model_set_parameters(fit0@Model, x)
    lavimplied <- lavaan::lav_model_implied(lavmodel)
    GLIST <- lavmodel@GLIST
    list2env(GLIST, environment())
    mu <- lavimplied$mean[[1]]
    Sigma <- lavimplied$cov[[1]]
  } else {
    # params in this order
    #
    # lambda: fx=~x2 fx=~x3 fy=~y2 fy=~y3
    # beta:   fy~fx
    # theta:  x1~~x1 x2~~x2 x3~~x3 y1~~y1 y2~~y2 y3~~y3
    # psi:    fx~~fx fy~~fy
    # nu:     x1~1   x2~1   x3~1   y1~1   y2~1   y3~1

    lambda <- matrix(c(
      1,    0,
      x[1], 0,
      x[2], 0,
      0,    1,
      0,    x[3],
      0,    x[4]
    ), ncol = 2, byrow = TRUE)
    beta <- matrix(c(0, x[5], 0, 0), nrow = 2, ncol = 2)
    theta <- diag(x[6:11])
    psi <- diag(x[12:13])

    alpha <- rep(0, 2)
    nu <- rep(0, 6)

    mu <- apply(data, 2, mean)
    Btil <- solve(diag(2) - beta)
    Sigma <- lambda %*% Btil %*% psi %*% t(Btil) %*% t(lambda) + theta
  }

  # Extra bits
  Sigma_inv <- solve(Sigma)
  Psitil <- Btil %*% psi %*% t(Btil)

  list(alpha = alpha, lambda = lambda, nu = nu, psi = psi, theta = theta,
       beta = beta, mu = mu, Sigma = Sigma,
       Sigma_inv = Sigma_inv, Btil = Btil, Psitil = Psitil)
}

LOGLIK_TWOFAC <- function(x, model, data) {
  with(make_matrices_twofac(x, data = data), {
    cen_dat <- scale(data, center = TRUE, scale = FALSE)
    sum(mvtnorm::dmvnorm(cen_dat, sigma = Sigma, log = TRUE))
  })
}

GRAD_TWOFAC <- function(x, model, data, i = NULL) {
  with(make_matrices_twofac(x, data = data), {

    Y <- as.matrix(data)
    if (!is.null(i)) Y <- matrix(Y[i, ], nrow = 1)
    n <- nrow(Y)

    # sample moments
    res <- sweep(Y, 2, as.vector(mu), "-")
    S <- crossprod(res) / n

    # variance components for gradient
    L <- Sigma_inv %*% (S - Sigma) %*% Sigma_inv
    lamLlam <- t(lambda) %*% L %*% lambda

    grad_theta <- (n / 2) * diag(L)
    grad_Psi <- (n / 2) * diag(t(Btil) %*% lamLlam %*% Btil)
    grad_beta <- n * (t(Btil) %*% lamLlam %*% psi)[1, 2]
    grad_lambda <- n * (L %*% lambda %*% Psitil)
    grad_lambda <- grad_lambda[cbind(c(2, 3, 5, 6), c(1, 1, 2, 2))]

    c(grad_lambda, grad_beta, grad_theta, grad_Psi)
  })
}

EMAT_TWOFAC <- function(x, model, data) {
  tcrossprod(sapply(seq_len(nrow(data)), \(i) GRAD_TWOFAC(x, model, data, i)))
}

JMAT_TWOFAC <- function(x, model, data) {
  -numDeriv::jacobian(GRAD_TWOFAC, x, model = model, data = data, i = NULL)
}

PENALTY_TWOFAC <- function(x, model, data, fallback = -10^8) {
  e <- EMAT_TWOFAC(x, model, data)
  jinv <- solve(JMAT_TWOFAC(x, model, data))
  -sum(jinv * e) / 2
}

BIAS_TWOFAC <- function(x, model, data) {
  jinv <- solve(JMAT_TWOFAC(x, model, data))
  A <- numDeriv::grad(PENALTY_TWOFAC, x, model = model, data = data)
  -drop(jinv %*% A)
}

#' Fit two factor SEM using lavaan syntax manually
#'
#' @inherit fit_sem params return
#' @param start A numeric vector of starting values
#'
#' @return A list with the following components:
#' \item{coefficients}{The estimated coefficients}
#' \item{stderr}{The standard errors of the estimated coefficients}
#' \item{bias}{The bias of the estimated coefficients}
#' \item{timing}{The elapsed time}
#' \item{converged}{Whether the optimization converged}
#' \item{optim}{The optimization object}
#'
#' @export
fit_twofac <- function(
    model,
    data,
    rbm = c("none", "explicit", "implicit"),
    # uselavaan = FALSE,
    start = NULL
) {
  rbm <- match.arg(rbm)
  if (is.null(start)) start <- c(1, 0, 1, 0, 0, 1)

  if (rbm == "none" | rbm == "explicit") {
    obj_fun <- function(x) -1 * LOGLIK_TWOFAC(x, model, data)
    grad_fun <- function(x) -1 * GRAD_TWOFAC(x, model, data)
  } else if (rbm == "implicit") {
    obj_fun <- function(x)
      -1 * (LOGLIK_TWOFAC(x, model, data) + PENALTY_TWOFAC(x, model, data))
    grad_fun <- NULL
  }

  start_time <- proc.time()

  res <- nlminb(start = start, objective = obj_fun, gradient = grad_fun)
  b <- if (rbm == "explicit") BIAS(res$par, model, data) else 0
  est <- res$par - b

  elapsed_time <- proc.time() - start_time
  elapsed_time <- elapsed_time["elapsed"]

  j <- JMAT_TWOFAC(est, model, data)
  jinv <- try(solve(j), silent = !TRUE)
  sds <- sqrt(diag(jinv))

  list(
    coefficients = est,
    stderr = sds,
    bias = b,
    timing = elapsed_time,
    converged = res$convergence == 0L,
    optim = res
  )

}
