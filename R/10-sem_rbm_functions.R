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
    plugin_pen = NULL,
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
  eigvals <- try(
    eigen(lavimplied$cov[[1]], symmetric = TRUE, only.values = TRUE)$values,
    silent = TRUE
  )
  # Sigma <- lavimplied$cov[[1]]  # THIS IS THE "current" Sigma

  if (any(eigvals < 1e-07) | inherits(eigvals, "try-error")) {
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

  if (!is.null(plugin_pen)) {
    if (!is.function(plugin_pen)) {
      cli::cli_abort("`plugin_pen` must be a function.")
    }

    pen_args <- list(
      x = theta,
      x.packed = theta.packed,
      lb = bounds$lb,
      ub = bounds$ub
    )
    pen_term <- do.call(plugin_pen, pen_args)
  }
  pen_term <- pen_term / lavsamplestats@ntotal

  # Verbose
  if (isTRUE(verbose)) {
    cat(
      "Log-lik =", loglik, "|",
      "Bias term =", bias_term, "|",
      "Penalty term =", pen_term, "|",
      "TOTAL =", loglik + bias_term - pen_term,
      "\n"
    )
  }

  # NOTE: bias_term is already negative
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
    kind = c("observed", "expected", "first.order")
  ) {

  kind <- match.arg(kind, c("observed", "expected", "first.order"))

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
  if (kind == "first.order") {
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

# Bias reduction penalty term
penalty <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions,
    kind
  ) {

  e <- information_matrix(
    theta = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    kind = "first.order"
  )
  j <- information_matrix(
    theta = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    kind = kind
  )
  if (TRUE) {#(lavmodel@eq.constraints) {
    jinv <- lav_model_information_augment_invert(
      lavmodel = lavmodel,
      information = j,
      inverted = TRUE,
      check.pd = FALSE,
      use.ginv = TRUE,
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
    kind_outside,
    kind_inside
  ) {

  j <- information_matrix(
    theta = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    kind = kind_outside
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
    lavoptions = lavoptions,
    kind = kind_inside
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
#' @inherit brlavaan params return
#' @param rbm The type of RBM method to use. One of `"none"`, `"explicit"`, or
#'   `"implicit"` (although, short forms are accepted, e.g. `"exp"`, `"iRBM`,
#'   etc.)
#' @param plugin_pen The type of penalty to use. One of `NULL`, `"pen_ridge"`,
#'   or `"pen_ridge_bound"`. User-specified penalty functions are possible. See
#'   description for details.
#' @param info_pen The type of information matrix to use for the penalty term.
#' @param info_bias The type of information matrix to use for the bias term of
#'   the explicit reduced bias method.
#' @param info_se The type of information matrix to use for calculation of
#'   standard errors.
#' @param debug If TRUE, the function will return a list of intermediate results
#'   for debugging purposes.
#' @param maxgrad If TRUE, the function will return the maximum gradient value.
#'
#' @export
fit_sem <- function(
    model,
    data,
    rbm = c("implicit", "explicit", "none"),
    plugin_pen = NULL,
    info_pen = "observed",
    info_bias = "observed",
    info_se = "observed",  #FIXME: this is akin to 'information' in lavaan
    debug = FALSE,
    lavfun = "sem",
    maxgrad = FALSE,
    ...
  ) {

  # Check arguments and choose method-------------------------------------------
  lavargs <- list(...)
  if ("method" %in% names(lavargs)) {
    lifecycle::deprecate_warn(
      when = "0.0.1",
      what = I("The `method` argument"),
      details = "Use `estimator.args` list to specify RBM options instead."
    )

    rbm <- lavargs$method
    plugin_pen <- NULL
    if (rbm == "ML") rbm <- "none"
    if (grepl("eRB|eBR", rbm, ignore.case = TRUE)) rbm <- "explicit"
    if (grepl("iRB|iBR", rbm, ignore.case = TRUE)) rbm <- "implicit"
    if (grepl("iRB|iBR", rbm, ignore.case = TRUE)) rbm <- "implicit"
    if (rbm == "iRBMp" | rbm == "iBRMp") {
      rbm <- "implicit"
      plugin_pen <- pen_huber
      cli::cli_alert_warning("Using the Huber penalty for iRBM.")
    }
    lavargs$method <- NULL

  } else if ("estimator.args" %in% names(lavargs)) {
    # This comes from the brlavaan() function, which overwrites the arguments
    if (!is.null(lavargs$estimator.args$rbm))
      rbm <- lavargs$estimator.args$rbm
    if (!is.null(lavargs$estimator.args$plugin_pen))
      plugin_pen <- lavargs$estimator.args$plugin_pen
    if (!is.null(lavargs$estimator.args$info_pen))
      info_pen <- lavargs$estimator.args$info_pen
    if (!is.null(lavargs$estimator.args$info_bias))
      info_bias <- lavargs$estimator.args$info_bias
  }
  if ("information" %in% names(lavargs)) info_se <- lavargs$information
  if ("estimator" %in% names(lavargs)) estimator <- lavargs$estimator
  else estimator <- "ML"

  # Validate arguments
  rbm <- validate_rbm(rbm)
  information <- list(info_pen, info_bias, info_se)
  names(information) <- c("penalty", "bias", "se")

  if (estimator != "ML") cli::cli_abort("Bias reduction methods are currently only available for ML estimation.")

  # Which method?
  is_ML     <- rbm == "none"
  is_eRBM   <- rbm == "explicit"
  is_iRBM   <- rbm == "implicit"

  # Initialise {lavaan} model object -------------------------------------------
  lavargs$model <- model
  lavargs$data <- data
  lavargs$do.fit <- FALSE

  fit0           <- do.call(get(lavfun, envir = asNamespace("lavaan")), lavargs)
  lavmodel       <- fit0@Model
  lavsamplestats <- fit0@SampleStats
  lavdata        <- fit0@Data
  lavoptions     <- fit0@Options
  n              <- lavaan::nobs(fit0)

  # Bounds NOTE: Since lavaan 0.6-19, settings bounds = TRUE will remove the
  # simple equality constraints. So, have to redo.
  lavargs$bounds <- TRUE  # bounds are always used
  fit0 <- do.call(get(lavfun, envir = asNamespace("lavaan")), lavargs)
  pt <- lavaan::partable(fit0)
  lb <- pt$lower[pt$free > 0]
  ub <- pt$upper[pt$free > 0]
  bounds <- list(lb = lb, ub = ub)

  # Starting values
  start <- lavargs$start
  if (is.null(start)) start <- lavaan::coef(fit0)  # starting values
  # Fix starting values, so they are not at boundaries
  fixidx <- which(start >= ub | start <= lb)
  start[fixidx] <- ((ub + lb) / 2)[fixidx]

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
      bounds = bounds,
      n = n,

      # Information about constraints
      lavmodel_eq.constraints = lavmodel@eq.constraints,
      lavmodel_eq.constraints.k0 = lavmodel@eq.constraints.k0,
      lavmodel_eq.constraints.K = lavmodel@eq.constraints.K,

      # Information about the model
      loglik = loglik(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions, bias_reduction = isFALSE(is_ML | is_eRBM), plugin_pen = plugin_pen, kind = info_pen, bounds = bounds, verbose = FALSE),
      grad_loglik = grad_loglik(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions),
      j = information_matrix(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions, kind = info_pen),
      jinv = lav_model_information_augment_invert(
        lavmodel = lavmodel,
        information = information_matrix(
          theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions, info_pen
        ),
        inverted = TRUE,
        check.pd = TRUE,
        use.ginv = FALSE,
        rm.idx = integer(0L)
      ),
      e = information_matrix(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions, kind = "first.order"),
      scores_loglik = scores_loglik(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions),
      penalty = penalty(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions, info_pen),
      bias = bias(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions, info_pen, info_bias),
      information = unlist(information)
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
      plugin_pen = plugin_pen,
      kind = info_pen,
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
      lavoptions = lavoptions,
      kind_inside = info_pen,
      kind_outside = info_bias
    )
  est <- res$par - b

  elapsed_time <- proc.time() - start_time
  elapsed_time <- elapsed_time["elapsed"]

  # Get max grad ---------------------------------------------------------------
  if (isTRUE(maxgrad)) {
    max_grad <- max(abs(numDeriv::grad(
      func = obj_fun,
      x = res$par,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )))
  } else {
    max_grad <- NULL
  }


  # Standard errors ------------------------------------------------------------
  j <- information_matrix(
    theta = as.numeric(est),
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    kind = info_se
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
    names(est) <- names(start)
  }

  list(
    coefficients = est,
    stderr = sds,
    timing = elapsed_time,
    converged = res$convergence == 0L,
    max_grad = max_grad,
    optim = res,
    vcov = jinv,
    information = information,
    lavfun = lavfun,
    estimator = estimator,
    rbm = rbm,
    plugin_pen = plugin_pen
  )
}
