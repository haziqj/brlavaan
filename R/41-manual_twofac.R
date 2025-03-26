# Adapting the manual script for the growth model by HJ to the factor model
# Manually obtaining derivatives for the two factor model to check correctness

make_matrices_twofac <- function(x, uselavaan = FALSE, model, data) {
  # model and data arguments used only for creating fit0
  if (isTRUE(uselavaan)) {
    # TODO: Fix meanstructure so using lavaan provides mu

    fit0 <- lavaan::sem(model = model, data = data, do.fit = FALSE, ceq.simple = TRUE)
    lavmodel <- lavaan::lav_model_set_parameters(fit0@Model, x)
    lavimplied <- lavaan::lav_model_implied(lavmodel)
    GLIST <- lavmodel@GLIST
    list2env(GLIST, environment())
    mu <- lavimplied$mean[[1]]
    Sigma <- lavimplied$cov[[1]]
  } else {
    nu <- colMeans(data)
    lambda <- matrix(c(
      1, 0,
      x[1], 0,
      x[2], 0,
      0, 1,
      0, x[3],
      0, x[4]
    ), ncol = 2, byrow = TRUE)
    beta <- matrix(c(0, x[5], 0, 0), nrow=2, ncol=2)
    psi <- matrix(c(x[1], x[5], x[5], x[3]), ncol = 2)
    theta <- diag(x[6:11])
    psi <- diag(x[12:13])
    IB_inv <- solve(diag(2) - beta)

    mu <- nu
    Sigma <- (lambda %*% IB_inv %*% psi %*% t(IB_inv) %*% t(lambda)) + theta
  }


  list(lambda = lambda, beta = beta, nu = nu, psi = psi, theta = theta,
       mu = mu, Sigma = Sigma)
}

LOGLIK_twofac <- function(x, model, data) {
  with(make_matrices_twofac(x, data=data), {
    sum(mvtnorm::dmvnorm(data, mean = mu, sigma = Sigma, log = TRUE))
  })
}

GRAD_twofac <- function(x, model, data) {
  with(make_matrices_twofac(x, data=data), {
    Sigma_inv <- solve(Sigma)
    IB_inv <- solve(diag(2) - beta)
    Y <- as.matrix(data)
    n <- nrow(Y)

    # sample moments
    res <- sweep(Y, 2, as.vector(mu), "-") # (Here mu is Ybar)
    S <- crossprod(res) / n

    # variance components for gradient
    E <- (n / 2) * Sigma_inv %*% (S - Sigma) %*% Sigma_inv
    M <- IB_inv %*% psi %*% t(IB_inv)

    grad_theta <- diag(E)
    grad_Psi <- diag(t(IB_inv) %*% t(lambda) %*% E %*% lambda %*% IB_inv)

    parts <- M %*% t(lambda) %*% (2*E)
    grad_lambda <-c(parts[1,2], parts[1,3], parts[2,5], parts[2,6])

    parts_b <- psi %*% (diag(2) + t(beta)) %*% t(lambda) %*% (2*E) %*% lambda
    grad_beta <- parts_b[1, 2]

    c(grad_lambda, grad_beta, grad_theta, grad_Psi)
  })
}






#set.seed(26)
#dat <- gen_data_twofac(n = 15, rel = 0.5, dist = "Normal")
#mod <- txt_mod_twofac(0.5)
#tru <- truth(dat)

#round(GRAD_twofac(tru, mod, dat) - numDeriv::grad(LOGLIK_twofac, tru, model=mod, data=dat), 5)
#GRAD_twofac(tru, mod, dat)
#numDeriv::grad(LOGLIK_twofac, tru, model=mod, data=dat)
