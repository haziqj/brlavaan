# Note: lavaan unexported functions are declared in 02-lavaan_unexported.R

# Log-likelihood function
loglik <- function(
    theta,  # packed / short version
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions,
    bias_reduction,
    kind,
    plugin_penalty = NULL,
    bounds,
    verbose
  ) {

  # Unpack theta
  theta.packed <- theta
  if (lavmodel@eq.constraints) {
    theta <- as.numeric(lavmodel@eq.constraints.K %*% theta) +
      lavmodel@eq.constraints.k0
  }

  # Fill in parameters in lavaan's internal matrix representation
  this_lavmodel <- lavaan::lav_model_set_parameters(lavmodel, x = theta)

  # Compute model implied mean and (co)variance matrix
  lavimplied <- lavaan::lav_model_implied(this_lavmodel)

  # Check if lavimplied$cov[[1]] is PD
  eigvals <- eigen(lavimplied$cov[[1]], symmetric = TRUE, only.values = TRUE)$values
  if (any(eigvals < 1e-07)) {
    # Return huge number, to signal the nlminb() optimizer something is not
    # quite OK with this parameter vector
    return(-1e40)
  }

  # Mean structure?
  if (lavmodel@meanstructure) {
    Mu <- lavimplied$mean[[1]]
  } else {
    Mu <- lavsamplestats@mean[[1]]
  }

  # Total log-likelihood value
  loglik <- lav_mvnorm_loglik_samplestats(
    sample.mean = lavsamplestats@mean[[1]],
    sample.cov  = lavsamplestats@cov[[1]],
    sample.nobs = lavsamplestats@nobs[[1]],
    Mu          = Mu,
    Sigma       = lavimplied$cov[[1]],
    x.idx       = lavsamplestats@x.idx[[1]],
    x.mean      = lavsamplestats@mean.x[[1]],
    x.cov       = lavsamplestats@cov.x[[1]],
    Sinv.method = "eigen",
    Sigma.inv   = NULL
  )

  # Bias reduction and plugin penalty term
  bias_term <- pen_term <- 0
  if (isTRUE(bias_reduction)) {
    bias_term <- penalty(
      theta = theta.packed,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      kind = kind
    )
  }

  if (!is.null(plugin_penalty)) {
    if (!is.function(plugin_penalty)) {
      cli::cli_abort("`plugin_penalty` must be a function.")
    }

    pen_args <- list(
      x = theta,
      x.packed = theta.packed,
      lb = bounds$lb,
      ub = bounds$ub
    )
    pen_term <- do.call(plugin_penalty, pen_args)
  }
  pen_term <- pen_term / lavsamplestats@ntotal

  # Verbose
  if (isTRUE(verbose)) {
    cat(
      "Log-lik =", loglik, "|",
      "Bias term =", bias_term, "|",
      "Penalty term =", pen_term, "|",
      "\n"
    )
  }

  out <- loglik + bias_term - pen_term
  # nlminb() works better with smaller numbers
  # out <- out / (2 * lavsamplestats@ntotal)

  out
}

# Gradient function
grad_loglik <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions
  ) {

  # Unpack theta
  if (lavmodel@eq.constraints) {
    theta <- as.numeric(lavmodel@eq.constraints.K %*% theta) +
      lavmodel@eq.constraints.k0
  }

  # Fill in parameters in lavaan's internal matrix representation
  this_lavmodel <- lavaan::lav_model_set_parameters(lavmodel, x = theta)

  # Gradient of fit function F_ML (not loglik yet)
  grad_F <- lav_model_gradient(
    lavmodel = this_lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata
  )
  # Rescale so we get gradient of loglik
  out <- -1 * lavsamplestats@ntotal * grad_F

  # Repack gradients
  if (lavmodel@eq.constraints) {
    out <- as.numeric(out %*% lavmodel@eq.constraints.K)
  }

  # Adjust gradient for rescaling
  # out <- out / (2 * lavsamplestats@ntotal)
  out
}

# Casewise scores (gradient of loglik for a single observation)
scores_loglik <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions
) {

  # Unpack theta
  if (lavmodel@eq.constraints) {
    theta <- as.numeric(lavmodel@eq.constraints.K %*% theta) +
      lavmodel@eq.constraints.k0
  }

  # Fill in parameters in lavaan's internal matrix representation
  this_lavmodel <- lavaan::lav_model_set_parameters(lavmodel, x = theta)

  # Compute model implied mean and (co)variance matrix
  lavimplied <- lavaan::lav_model_implied(this_lavmodel)

  # Only 1 group
  moments <- list(cov = lavimplied$cov[[1]])
  if(lavmodel@meanstructure) {
    moments$mean <- lavimplied$mean[[1]]
  }

  ntab <- unlist(lavdata@norig)
  ntot <- sum(ntab)
  npar <- length(theta)
  out <- lav_scores_ml(
    ntab = ntab,
    ntot = ntot,
    npar = npar,
    moments = moments,
    lavdata = lavdata,
    lavsamplestats = lavsamplestats,
    lavmodel = lavmodel,
    lavoptions = lavoptions,
    scaling = FALSE
  )

  out
}

# Hessian function
hessian_loglik <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions
  ) {

  # Unpack theta
  if (lavmodel@eq.constraints) {
    theta <- as.numeric(lavmodel@eq.constraints.K %*% theta) +
      lavmodel@eq.constraints.k0
  }

  # Fill in parameters in lavaan's internal matrix representation
  this_lavmodel <- lavaan::lav_model_set_parameters(lavmodel, x = theta)

  # Hessian of F_ML
  out <- lav_model_hessian(
    lavmodel = this_lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )

  # Rescale so we get Hessian of loglik
  -1 * lavsamplestats@ntotal * out
}

# Full information matrix
information_matrix <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions,
    kind = c("observed", "expected", "firstorder")
  ) {

  kind <- match.arg(kind, c("observed", "expected", "firstorder"))

  # Unpack theta
  if (lavmodel@eq.constraints) {
    theta <- as.numeric(lavmodel@eq.constraints.K %*% theta) +
      lavmodel@eq.constraints.k0
  }

  # Fill in parameters in lavaan's internal matrix representation
  this_lavmodel <- lavaan::lav_model_set_parameters(lavmodel, x = theta)

  # Get information matrix
  if (kind == "observed") {
    # FIXME: When is this needed?
    # Change to 'simplified' Hessian
    if (FALSE) lavoptions$observed.information <- c("h1", "h1")
    out <- lav_model_information_observed(
      lavmodel = this_lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
  }
  if (kind == "expected") {
    out <- lav_model_information_expected(
      lavmodel = this_lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
  }
  if (kind == "firstorder") {
    out <- lav_model_information_firstorder(
      lavmodel = this_lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
  }

  # This is unit information, so multiply by N
  lavsamplestats@ntotal * out
}

# Observed information
information_observed <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions
  ) {
  information_matrix(
    theta = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    kind = "observed"
  )
}

# Expected information
information_expected <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions
  ) {
  information_matrix(
    theta = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    kind = "expected"
  )
}

# Outer product of casewise scores
information_firstorder <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions
  ) {
  information_matrix(
    theta = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    kind = "firstorder"
  )
}

# Bias reduction penalty term
penalty <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions,
    kind = c("observed", "expected")
  ) {

  kind <- match.arg(kind, c("observed", "expected"))

  e <- information_firstorder(
    theta = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
  j <- information_matrix(
    theta = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    kind = kind
  )
  if (lavmodel@eq.constraints) {
    jinv <- lav_model_information_augment_invert(
      lavmodel = lavmodel,
      information = j,
      inverted = TRUE,
      check.pd = FALSE,
      use.ginv = FALSE,
      rm.idx = integer(0L)
    )
  } else {
    jinv <- solve(j)
  }

  if (inherits(jinv, "try-error")) {
    return(NA)
  } else {
    return(-sum(jinv * e) / 2)
  }
}

# The bias term for explicit bias reduction method
bias <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions,
    kind = c("observed", "expected")
  ) {

  kind <- match.arg(kind, c("observed", "expected"))

  j <- information_matrix(
    theta = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    kind = kind
  )
  jinv <- lav_model_information_augment_invert(
    lavmodel = lavmodel,
    information = j,
    inverted = TRUE,
    check.pd = FALSE,
    use.ginv = FALSE,
    rm.idx = integer(0L)
  )
  A <- numDeriv::grad(
    func = penalty,
    x = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )

  # Unpack A
  if (lavmodel@eq.constraints) {
    A <- lavmodel@eq.constraints.K %*% A
  }

  out <- -drop(jinv %*% A)

  # Repack! (since optimiser expects packed version)
  if (lavmodel@eq.constraints) {
    out <- as.numeric(out %*% lavmodel@eq.constraints.K)
  }

  out
}

#' Fit a structural equation model using bias-reducing methods
#'
#' This is the workhorse function to fit a structural equation model using
#' bias-reducing methods which we use for our simulations. For endusers, we
#' recommend using the [brsem()], [brcfa()], or [brgrowth()] functions which is
#' a wrapper around this function and outputs something more user-friendly.
#'
#' @inherit brsem params return
#' @param debug If TRUE, the function will return a list of intermediate
#'   results for debugging purposes.
#' @param lavfun The lavaan function to use. Default is "sem".
#'
#' @export
fit_sem <- function(
    model,
    data,
    estimator = "ML",
    estimator.args = list(rbm = FALSE, plugin_penalty = NULL),
    information = c("expected", "observed", "first.order"),
    debug = FALSE,
    lavfun = "sem",
    ...
  ) {

  information <- match.arg(information, c("expected", "observed", "first.order"))
  orig_info <- information
  if (orig_info == "first.order") {
    information <- "expected"
    cli::cli_alert_info("Bias reduction methods will use expected information matrix, and standard error computation will use the outer product of the casewise scores.")
  }

  # Initialise {lavaan} model object -------------------------------------------
  lavargs <- list(...)
  lavargs$model <- model
  lavargs$data <- data
  lavargs$do.fit <- FALSE

  if (estimator != "ML") {
    cli::cli_abort("Bias reduction methods are currently only available for ML estimation.")
  }

  # Catch old arguments
  if ("method" %in% names(lavargs)) {
    # cli::cli_alert_warning("RBM options are now specified in 'estimator.args' list  instead of 'method' argument.")
    lifecycle::deprecate_warn(
      when = "0.0.1",
      what = I("The `method` argument"),
      details = "Use `estimator.args` list to specify RBM options instead."
    )

    rbm <- lavargs$method
    plugin_penalty <- NULL
    if (rbm == "ML") rbm <- FALSE
    if (rbm == "iRBMp") {
      rbm <- "iRBM"
      plugin_penalty <- pen_huber
      cli::cli_alert_warning("Using the Huber penalty for iRBM.")
    }
    lavargs$method <- NULL
  } else {
    rbm <- estimator.args$rbm
    plugin_penalty <- estimator.args$plugin_penalty
  }
  if (!isFALSE(rbm)) rbm <- match.arg(rbm, c("eRBM", "iRBM"))

  # Which method?
  is_ML     <- isFALSE(rbm)
  is_eRBM   <- rbm == "eRBM"
  is_iRBM   <- rbm == "iRBM"

  fit0           <- do.call(get(lavfun, envir = asNamespace("lavaan")), lavargs)
  lavmodel       <- fit0@Model
  lavsamplestats <- fit0@SampleStats
  lavdata        <- fit0@Data
  lavoptions     <- fit0@Options

  start <- lavargs$start
  if (is.null(start)) start <- lavaan::coef(fit0)  # starting values
  n <- lavaan::nobs(fit0)

  # Bounds NOTE: Since lavaan 0.6-19, settings bounds = TRUE will remove the
  # simple equality constraints. So, have to redo.
  lavargs$bounds <- TRUE  # bounds are always used
  fit0           <- do.call(get(lavfun, envir = asNamespace("lavaan")), lavargs)
  pt <- lavaan::partable(fit0)
  lb <- pt$lower[pt$free > 0]
  ub <- pt$upper[pt$free > 0]
  bounds <- list(lb = lb, ub = ub)

  # Pack theta and bounds
  if (lavmodel@eq.constraints) {
    theta_pack <- as.numeric(
        (start - lavmodel@eq.constraints.k0) %*% lavmodel@eq.constraints.K
      )
  } else {
    theta_pack <- start
  }

  # DEBUG ----------------------------------------------------------------------
  if (isTRUE(debug)) {
    out <- list(
      # Parameters and sample size
      theta_unpack = start,
      theta_pack = theta_pack,
      n = n,

      # Information about constraints
      lavmodel_eq.constraints = lavmodel@eq.constraints,
      lavmodel_eq.constraints.k0 = lavmodel@eq.constraints.k0,
      lavmodel_eq.constraints.K = lavmodel@eq.constraints.K,

      # Information about the model
      loglik = loglik(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions, bias_reduction = isFALSE(is_ML | is_eRBM), plugin_penalty = plugin_penalty, kind = information, bounds = bounds, verbose = FALSE),
      grad_loglik = grad_loglik(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions),
      j = information_matrix(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions, kind = orig_info),
      jinv = lav_model_information_augment_invert(
        lavmodel = lavmodel,
        information = information_matrix(
          theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions, orig_info
        ),
        inverted = TRUE,
        check.pd = FALSE,
        use.ginv = FALSE,
        rm.idx = integer(0L)
      ),
      e = information_firstorder(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions),
      scores_loglik = scores_loglik(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions),
      penalty = penalty(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions, information),
      bias = bias(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions, information),
      information_penalty = information,
      information_se = orig_info
    )
    return(out)
  }

  # Verbose
  verbose <- lavargs$verbose
  if (is.null(verbose)) verbose <- FALSE

  # Start fitting process ------------------------------------------------------
  start_time <- proc.time()
  obj_fun <- function(x, ...) {
    -1 * loglik(
      theta = x,
      ...,
      bias_reduction = isFALSE(is_ML | is_eRBM),
      plugin_penalty = plugin_penalty,
      kind = information,
      bounds = bounds,
      verbose = verbose
    )
  }
  grad_fun <-
    if (isTRUE(is_ML | is_eRBM))
      function(x, ...) -1 * grad_loglik(x, ...)
    else
      NULL

  res <- nlminb(
    start = theta_pack,
    objective = obj_fun,
    gradient = grad_fun,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    control = lavargs$control
  )

  b <- 0
  if (is_eRBM)
    b <- b + bias(
      theta = res$par,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
  est <- res$par - b

  elapsed_time <- proc.time() - start_time
  elapsed_time <- elapsed_time["elapsed"]

  # Standard errors ------------------------------------------------------------
  j <- information_matrix(
    theta = as.numeric(est),
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    kind = orig_info
  )
  jinv <- lav_model_information_augment_invert(
    lavmodel = lavmodel,
    information = j,
    inverted = TRUE,
    check.pd = FALSE,
    use.ginv = FALSE,
    rm.idx = integer(0L)
  )
  sds <- sqrt(diag(jinv))

  # Unpack estimators
  if (lavmodel@eq.constraints) {
    est  <- as.numeric(lavmodel@eq.constraints.K %*% est) +
      lavmodel@eq.constraints.k0
  }

  list(
    coefficients = est,
    stderr = sds,
    timing = elapsed_time,
    converged = res$convergence == 0L,
    optim = res,
    vcov = jinv,
    information_penalty = information,
    information_se = orig_info,
    lavfun = lavfun,
    estimator = estimator,
    rbm = rbm,
    plugin_penalty = plugin_penalty
  )
}
