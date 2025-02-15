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
    verbose,
    nearPD = TRUE,
    fn.scale = 1
  ) {

  # Unpack theta
  theta.packed <- theta
  theta <- unpack_x(theta.packed, lavmodel)

  # Fill in parameters in lavaan's internal matrix representation
  this_lavmodel <- lavaan::lav_model_set_parameters(lavmodel, x = theta)

  # Compute model implied mean and (co)variance matrix
  lavimplied <- lavaan::lav_model_implied(this_lavmodel)

  # Check if implied Sigma is PD
  Sigma <- lavimplied$cov[[1]]  # THIS IS THE "current" Sigma
  eigvals <- try(
    eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values,
    silent = TRUE
  )

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
    Sigma       = Sigma,
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
      kind = kind,
      nearPD = nearPD
    )
  }

  if (!is.null(plugin_pen)) {
    if (!is.function(plugin_pen)) {
      cli::cli_abort("`plugin_pen` must be a function.")
    }

    pen_args <- list(
      x = theta,
      x.packed = theta.packed,
      lb = bounds$lower,
      ub = bounds$upper
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
  out <- fn.scale * out

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
  theta <- unpack_x(theta, lavmodel)

  # Update lavmodel
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
  out <- pack_x(out, lavmodel)

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
  theta <- unpack_x(theta, lavmodel)

  # Update lavstuff
  this_lavmodel <- lavaan::lav_model_set_parameters(lavmodel, x = theta)
  lavimplied <- lavaan::lav_model_implied(this_lavmodel)

  # FIXME: Only 1 group!!!
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
  theta <- unpack_x(theta, lavmodel)

  # Update lavstuff
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
  theta <- unpack_x(theta, lavmodel)

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
    kind,
    nearPD = TRUE
  ) {

  e <- information_matrix(
    theta = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    kind = "first.order"
  )
  if (isTRUE(nearPD)) e <- as.matrix(Matrix::nearPD(e)$mat)
  j <- information_matrix(
    theta = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    kind = kind
  )
  if (isTRUE(nearPD)) j <- as.matrix(Matrix::nearPD(j)$mat)
  jinv <- try(solve(j), silent = !TRUE)

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
    kind_inside,
    nearPD = TRUE
  ) {

  j <- information_matrix(
    theta = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    kind = kind_outside
  )
  if (isTRUE(nearPD)) j <- as.matrix(Matrix::nearPD(j)$mat)
  jinv <- try(solve(j), silent = !TRUE)

  A <- numDeriv::grad(
    func = penalty,
    x = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    kind = kind_inside,
    nearPD = nearPD
  )

  A <- unpack_x(A, lavmodel)  # unpack A
  out <- -drop(jinv %*% A)
  out <- pack_x(out, lavmodel)  # repack! since optimiser expects packed version

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
#' @param nearPD `r lifecycle::badge("experimental")` If TRUE, the function will
#'   ensure the information matrix is positive definite (using
#'   [Matrix::nearPD()]), which may resolve numerical issues in the bias
#'   reduction method and standard error calculations.
#' @param fn.scale `r lifecycle::badge("experimental")` A scaling factor for the
#'   log-likelihood function.
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
    nearPD = TRUE,
    fn.scale = 1,
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
  if (!is.null(lavargs$bounds))
    lavargs$ceq.simple <- TRUE  # YR advises to force simple bounds
  lavargs$do.fit <- FALSE

  fit0           <- do.call(get(lavfun, envir = asNamespace("lavaan")), lavargs)
  lavmodel       <- fit0@Model
  lavsamplestats <- fit0@SampleStats
  lavdata        <- fit0@Data
  lavoptions     <- fit0@Options
  lavpartable    <- fit0@ParTable
  n              <- fit0@SampleStats@ntotal

  # Starting values
  start <- lavaan::lav_model_get_parameters(lavmodel)
  if (lavmodel@eq.constraints) start <- pack_x(start, lavmodel)

  # Bounds (taken from lavaan's lav_model_estimate.R lines 240-277)
  if (is.null(lavpartable$lower)) {
    lower <- -Inf
  } else {
    if (lavmodel@ceq.simple.only) {
      free.idx <- which(lavpartable$free > 0L & !duplicated(lavpartable$free))
      lower <- lavpartable$lower[free.idx]
    } else if (lavmodel@eq.constraints) {
      # Bounds have no effect any longer...
      # 0.6-19 -> we switch to (un)constrained estimation
      lower <- -Inf
    } else {
      lower <- lavpartable$lower[lavpartable$free > 0L]
    }
  }
  if (is.null(lavpartable$upper)) {
    upper <- +Inf
  } else {
    if (lavmodel@ceq.simple.only) {
      free.idx <- which(lavpartable$free > 0L & !duplicated(lavpartable$free))
      upper <- lavpartable$upper[free.idx]
    } else if (lavmodel@eq.constraints) {
      # Bounds have no effect any longer...
      upper <- +Inf
    } else {
      upper <- lavpartable$upper[lavpartable$free > 0L]
    }
  }

  # DEBUG ----------------------------------------------------------------------
  if (isTRUE(debug)) {
    out <- list(
      # Parameters and sample size
      theta_unpack = unpack_x(start, lavmodel),
      theta_pack = start,
      lower = lower,
      upper = upper,
      n = n,

      # Information about constraints
      lavmodel_eq.constraints = lavmodel@eq.constraints,
      lavmodel_eq.constraints.k0 = lavmodel@eq.constraints.k0,
      lavmodel_eq.constraints.K = lavmodel@eq.constraints.K,

      # Information about the model
      loglik = loglik(start, lavmodel, lavsamplestats, lavdata, lavoptions,
                      bias_reduction = isFALSE(is_ML | is_eRBM),
                      plugin_pen = plugin_pen, kind = info_pen,
                      bounds = list(lower, upper), verbose = FALSE,
                      fn.scale = fn.scale, nearPD = nearPD),
      grad_loglik = grad_loglik(start, lavmodel, lavsamplestats, lavdata,
                                lavoptions),
      j = information_matrix(start, lavmodel, lavsamplestats, lavdata,
                             lavoptions, kind = info_pen),
      jinv = try(solve(
        information_matrix(start, lavmodel, lavsamplestats, lavdata, lavoptions,
                           kind = info_pen)
      )),
      e = information_matrix(start, lavmodel, lavsamplestats, lavdata,
                             lavoptions, kind = "first.order"),
      scores_loglik = scores_loglik(start, lavmodel, lavsamplestats, lavdata,
                                    lavoptions),
      penalty = penalty(start, lavmodel, lavsamplestats, lavdata, lavoptions,
                        info_pen, nearPD),
      bias = bias(start, lavmodel, lavsamplestats, lavdata, lavoptions,
                  info_pen, info_bias, nearPD),
      information = unlist(information),
      nearPD = nearPD
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
      bounds = list(lower = lower, upper = upper),
      verbose = verbose,
      nearPD = nearPD,
      fn.scale = fn.scale
    )
  }
  grad_fun <-
    if (isTRUE(is_ML | is_eRBM))
      function(x, ...) -1 * grad_loglik(x, ...)
    else
      NULL

  res <- nlminb(
    start = start,
    objective = obj_fun,
    gradient = grad_fun,
    lower = lower,
    upper = upper,
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
      kind_outside = info_bias,
      nearPD = nearPD
    )
  est <- res$par - b

  elapsed_time <- proc.time() - start_time
  elapsed_time <- elapsed_time["elapsed"]

  # Implied model covariance matrix --------------------------------------------
  lavmodel <- lavaan::lav_model_set_parameters(lavmodel, x = est)
  Sigma <- lavaan::lav_model_implied(lavmodel)$cov[[1]]

  # Get max grad ---------------------------------------------------------------
  if (isTRUE(maxgrad)) {
    hinv <- solve(numDeriv::hessian(
      func = obj_fun,
      x = res$par,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    ))
    max_scores <- numDeriv::grad(
      func = obj_fun,
      x = res$par,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
    scaled_grad <- as.numeric(hinv %*% max_scores)
  } else {
    scaled_grad <- NULL
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
  if (isTRUE(nearPD)) j <- as.matrix(Matrix::nearPD(j)$mat)
  jinv <- try(solve(j), silent = !TRUE)
  # if (lavmodel@eq.constraints) {
  #   jinv <- lav_model_information_augment_invert(
  #     lavmodel = lavmodel,
  #     information = j,
  #     inverted = TRUE,
  #     check.pd = FALSE,
  #     use.ginv = TRUE,
  #     rm.idx = integer(0L)
  #   )
  # }
  sds <- sqrt(diag(jinv))

  # Unpack estimators
  est <- unpack_x(est, lavmodel)
  names(est) <- names(lavaan::coef(fit0))

  list(
    coefficients = est,
    stderr = sds,
    timing = elapsed_time,
    converged = res$convergence == 0L,
    scaled_grad = scaled_grad,
    optim = res,
    vcov = jinv,
    Sigma = Sigma,
    information = information,
    lavfun = lavfun,
    estimator = estimator,
    rbm = rbm,
    plugin_pen = plugin_pen,
    bounds = list(lower = lower, upper = upper)
  )
}

Sigma_twofac <- function(x) {

  mod <- txt_mod_twofac(0.8)
  dat <- gen_data_twofac(100, 0.8)
  fit0 <- lavaan::sem(mod, dat)

  this_lavmodel <- lavaan::lav_model_set_parameters(fit0@Model, x = x)

  # Compute model implied mean and (co)variance matrix
  lavimplied <- lavaan::lav_model_implied(this_lavmodel)
  lavimplied$cov[[1]]
}

Sigma_growth <- function(x) {

  mod <- txt_mod_growth(0.8)
  dat <- gen_data_growth(100, 0.8)
  fit0 <- lavaan::sem(mod, dat)

  this_lavmodel <- lavaan::lav_model_set_parameters(fit0@Model, x = x)

  # Compute model implied mean and (co)variance matrix
  lavimplied <- lavaan::lav_model_implied(this_lavmodel)
  lavimplied$cov[[1]]
}

trf_x <- function(x, lavmodel, how = c("pack", "unpack")) {
  how <- match.arg(how)
  if (is.null(x)) return(x)
  if (lavmodel@eq.constraints) {
    k0 <- lavmodel@eq.constraints.k0
    K <- lavmodel@eq.constraints.K
    if (how == "pack") {
      out <- as.numeric((x - k0) %*% K)
    } else if (how == "unpack") {
      out <- as.numeric(K %*% x + k0)
    }
  } else if (lavmodel@ceq.simple.only) {
    K <- lavmodel@ceq.simple.K
    if (how == "pack") {
      out <- as.numeric(x %*% apply(K, 2, \(x) x / sum(x)))
    } else if (how == "unpack") {
      out <- as.numeric(K %*% x)
    }
  } else {
    out <- x
  }
  out
}

pack_x <- function(x, lavmodel) trf_x(x, lavmodel, "pack")
unpack_x <- function(x, lavmodel) trf_x(x, lavmodel, "unpack")
