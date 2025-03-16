make_matrices <- function(x, uselavaan = FALSE, model, data) {
  if (isTRUE(uselavaan)) {
    fit0 <- lavaan::growth(model = model, data = data, do.fit = FALSE, ceq.simple = TRUE)
    lavmodel <- lavaan::lav_model_set_parameters(fit0@Model, x)
    lavimplied <- lavaan::lav_model_implied(lavmodel)
    GLIST <- lavmodel@GLIST
    list2env(GLIST, environment())  # alpha, nu, lambda, psi, theta
    mu <- lavimplied$mean[[1]]
    Sigma <- lavimplied$cov[[1]]
  } else {
    alpha <- matrix(c(x[2], x[4]), ncol = 1)
    lambda <- matrix(c(
      1, 0,
      1, 1,
      1, 2,
      1, 3,
      1, 4,
      1, 5,
      1, 6,
      1, 7,
      1, 8,
      1, 9
    ), ncol = 2, byrow = TRUE)
    nu <- matrix(0, ncol = 1, nrow = nrow(lambda))
    psi <- matrix(c(x[1], x[5], x[5], x[3]), ncol = 2)
    theta <- diag(x[6], nrow = nrow(lambda))
    mu <- nu + lambda %*% alpha
    Sigma <- lambda %*% psi %*% t(lambda) + theta
  }

  list(alpha = alpha, lambda = lambda, nu = nu, psi = psi, theta = theta, mu = mu, Sigma = Sigma)
}

LOGLIK <- function(x, model = mod, data = dat) {
  list2env(make_matrices(x), environment())
  sum(mvtnorm::dmvnorm(dat, mean = mu, sigma = Sigma, log = TRUE))
}

GRAD <- function(x, model = mod, data = dat) {
  list2env(make_matrices(x), environment())
  Sigma_inv <- solve(Sigma)

  Y <- as.matrix(dat)
  n <- nrow(Y)
  p <- nrow(lambda)
  q <- ncol(lambda)

  # sample moments
  Ybar <- colMeans(Y)
  res <- sweep(Y, 2, as.vector(mu), "-")
  S <- crossprod(res) / n

  # mean components for gradient
  grad_mu <- Sigma_inv %*% (n * (Ybar - mu) / 2)
  grad_alpha <- t(lambda) %*% grad_mu

  # variance components for gradient
  E <- (n / 2) * Sigma_inv %*% (S - Sigma) %*% Sigma_inv
  grad_Psi <- t(lambda) %*% E %*% lambda
  grad_theta <- sum(diag(E))

  c(grad_Psi[1, 1], grad_alpha[1], grad_Psi[2, 2], grad_alpha[2], grad_Psi[1, 2], grad_theta)
}

EMAT <- function(x, model = mod, data = dat) {
  tcrossprod(sapply(seq_len(nrow(dat)), \(i) GRAD(x, model, data)))
}

JMAT <- function(x, model = mod, data = dat) {
  -numDeriv::hessian(LOGLIK, x, model = model, data = data)
}

PENALTY <- function(x, model = mod, data = dat) {
  e <- EMAT(x, model, data)
  jinv <- try(solve(JMAT(x, model, data)))
  if (inherits(jinv, "try-error")) {
    return(NA)
  } else {
    return(-sum(jinv * e) / 2)
  }
}

BIAS <- function(x, model = mod, data = dat) {
  jinv <- solve(JMAT(x, model, data))
  A <- numDeriv::grad(PENALTY, x, model = model, data = data)
  -drop(jinv %*% A)
}


#' Fit growth curve models using lavaan syntax manually
#'
#' @inherit fit_sem params return
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
fit_growth <- function(
    model,
    data,
    rbm = c("none", "explicit", "implicit"),
    # uselavaan = FALSE,
    start = NULL
  ) {
  rbm <- match.arg(rbm)
  if (is.null(start)) start <- c(1, 0, 1, 0, 0, 1)

  if (rbm == "none" | rbm == "explicit") {
    obj_fun <- function(x) -1 * LOGLIK(x, model, data)
    grad_fun <- function(x) -1 * GRAD(x, model, data)
  } else if (rbm == "implicit") {
    obj_fun <- function(x) {
      LOGLIK(x, model = model, data = data) +
        PENALTY(x, model = model, data = data)
    }
    grad_fun <- NULL
  }

  start_time <- proc.time()

  res <- nlminb(start = start, objective = obj_fun, gradient = grad_fun)
  b <- if (rbm == "explicit") BIAS(res$par, model, data) else 0
  est <- res$par - b

  elapsed_time <- proc.time() - start_time
  elapsed_time <- elapsed_time["elapsed"]

  j <- JMAT(est, model, data)
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
