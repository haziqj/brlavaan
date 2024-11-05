# Empirical bias correction ----------------------------------------------------
# 1. JJJ(theta) function is the variability matrix at theta
# 2. HHH(theta) function is the sensitivity matrix at theta
# 3. AAA(theta) function is the gradient of -1/2 tr(H^{-1} J}) at theta

JJJ <- function(theta, object = fit) {
  # extract slots
  lavmodel       <- object@Model
  lavsamplestats <- object@SampleStats
  lavdata        <- object@Data
  lavoptions     <- object@Options

  # update lavmodel with 'new' set of values
  lavmodel <- lav_model_set_parameters(lavmodel, x = theta)
  # FIXME: Is this correct? These would be the free parameters in the order of
  # partable(object)

  lavaan:::lav_model_information_firstorder(
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
}

HHH <- function(theta, object = fit, unit_information = FALSE) {
  # extract slots
  lavmodel       <- object@Model
  lavsamplestats <- object@SampleStats
  lavdata        <- object@Data
  lavoptions     <- object@Options
  lavcache       <- object@Cache
  lavimplied     <- object@implied
  lavh1          <- object@h1

  # update lavmodel with 'new' set of values
  lavmodel <- lav_model_set_parameters(lavmodel, x = theta)

  out <- lavaan:::lav_model_hessian(
    lavmodel       = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata        = lavdata,
    lavoptions     = lavoptions,
    lavcache       = lavcache,
    ceq.simple     = FALSE
  )

  if (isTRUE(unit_information)) out <- out / lavsamplestats@ntotal
  out
}

# Check
# J <- JJJ(theta = coef(fit))
# H <- HHH(theta = coef(fit))
# sqrt(diag(solve(H) %*% J %*% solve(H)) / n)  # sandwich se
# partable(fit)$se[partable(fit)$free > 0]

AAA <- function(theta, object = fit) {
  tmp <- function(x) {
    Hinv <- solve(HHH(theta = x, object = object, unit_information = FALSE))
    J    <- JJJ(x, object = object)
    -0.5 * sum(diag(Hinv %*% J))
  }
  numDeriv::grad(func = tmp, x = theta)
}

rb_empr <- function(fit) {
  start_time <- Sys.time()
  theta_hat <- coef(fit)
  n <- fit@Data@nobs[[1]]

  A <- AAA(theta_hat, fit)
  Hinv <- solve(HHH(theta_hat, fit, unit_information = FALSE))
  theta_tilde <- c(theta_hat + Hinv %*% A / n)  # FIXME: I think???

  end_time <- Sys.time()

  names(theta_tilde) <- names(coef(fit))
  attr(theta_tilde, "timing") <- end_time - start_time
  theta_tilde
}
