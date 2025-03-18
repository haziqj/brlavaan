# Note: lavaan unexported functions are declared in 02-lavaan_unexported.R

loglik <- function(
    x,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions,
    bias_reduction,
    plugin_pen = NULL,
    bounds,
    verbose,
    fn.scale = 1
  ) {

  lavmodel_x <- lavaan::lav_model_set_parameters(lavmodel, x)
  lavimplied <- lavaan::lav_model_implied(lavmodel_x)

  # Check if implied Sigma is PD
  Sigma <- lavimplied$cov[[1]]
  if (check_mat(Sigma)) return(-1e40)

  # Mean structure?
  if (lavmodel@meanstructure) {
    Mu <- lavimplied$mean[[1]]
  } else {
    Mu <- lavsamplestats@mean[[1]]
  }

  # Total log-likelihood value
  loglik <- lavaan___lav_mvnorm_loglik_samplestats(
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
      x = x,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
  }
  if (!is.null(plugin_pen)) {
    if (!is.function(plugin_pen)) {
      cli::cli_abort("`plugin_pen` must be a function.")
    }
    pen_args <- list(
      x = x,
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

  out <- loglik + bias_term - pen_term
  # nlminb() works better with smaller numbers
  # out <- out / (2 * lavsamplestats@ntotal)
  out <- fn.scale * out
  out
}

grad_loglik <- function(
    x,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions
  ) {

  # Gradient of fit function F_ML (not loglik yet)
  grad_F <- lavaan___lav_model_gradient(
    lavmodel = lavaan::lav_model_set_parameters(lavmodel, x),
    lavsamplestats = lavsamplestats,
    lavdata = lavdata
  )
  # Rescale so we get gradient of loglik
  out <- -1 * lavsamplestats@ntotal * grad_F

  # Repack gradients if needed
  if (lavmodel@ceq.simple.only)
    as.numeric(out %*% lavmodel@ceq.simple.K)
  else
    out
}

information_matrix <- function(
    x,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions,
    kind = c("observed", "expected", "first.order")
  ) {

  kind <- match.arg(kind)
  lavmodel_x <- lavaan::lav_model_set_parameters(lavmodel, x)
  lavimplied <- lavaan::lav_model_implied(lavmodel_x)

  if (kind == "observed") {
    # FIXME: When is this needed?
    # Change to 'simplified' Hessian
    if (FALSE) lavoptions$observed.information <- c("h1", "h1")
    out <- lavaan___lav_model_information_observed(
      lavmodel = lavmodel_x,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
  }
  if (kind == "expected") {
    out <- lavaan___lav_model_information_expected(
      lavmodel = lavmodel_x,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
  }
  if (kind == "first.order") {
    out <- lavaan___lav_model_information_firstorder(
      lavmodel = lavmodel_x,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
  }

  # This is unit information, so multiply by N
  lavsamplestats@ntotal * out
}

penalty <- function(
    x,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions,
    fallback = -10^8
  ) {

  e <- information_matrix(
    x = x,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    kind = "first.order"
  )
  j <- information_matrix(
    x = x,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    kind = "observed"
  )

  if (lavmodel@ceq.simple.only) {
    K <- lavmodel@ceq.simple.K
    e <- t(K) %*% e %*% K
    j <- t(K) %*% j %*% K
  }
  if (check_mat(j)) return(fallback)
  jinv <- try(solve(j), silent = TRUE)
  if (inherits(jinv, "try-error")) fallback else -sum(jinv * e) / 2

}

# The bias term for explicit bias reduction method
bias <- function(
    x,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions
  ) {

  j <- information_matrix(
    x = x,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    kind = "observed"
  )
  if (lavmodel@ceq.simple.only) {
    K <- lavmodel@ceq.simple.K
    j <- t(K) %*% j %*% K
  }
  if (check_mat(j)) return(rep(NA, length(x)))
  else jinv <- try(solve(j), silent = TRUE)

  A <- numDeriv::grad(
    func = penalty,
    x = x,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )

  -drop(jinv %*% A)
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
#' @param debug If TRUE, the function will return a list of intermediate results
#'   for debugging purposes.
#' @param maxgrad If TRUE, the function will return the maximum gradient value.
#' @param fn.scale `r lifecycle::badge("experimental")` A scaling factor for the
#'   log-likelihood function.
#' @export
fit_sem <- function(
    model,
    data,
    rbm = c("implicit", "explicit", "none"),
    plugin_pen = NULL,
    debug = FALSE,
    lavfun = "sem",
    maxgrad = FALSE,
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
  }
  if ("information" %in% names(lavargs)) information <- lavargs$information
  else information <- "observed"
  if ("estimator" %in% names(lavargs)) estimator <- lavargs$estimator
  else estimator <- "ML"

  # Validate arguments
  rbm <- validate_rbm(rbm)
  if (estimator != "ML") cli::cli_abort("Bias reduction methods are currently only available for ML estimation.")

  # Which method?
  is_ML     <- rbm == "none"
  is_eRBM   <- rbm == "explicit"
  is_iRBM   <- rbm == "implicit"

  # Initialise {lavaan} model object -------------------------------------------
  lavargs$model <- model
  lavargs$data <- data
  lavargs$ceq.simple <- TRUE  # force ceq.simple rather than eq.constraints
  lavargs$do.fit <- FALSE

  fit0           <- do.call(get(lavfun, envir = asNamespace("lavaan")), lavargs)
  lavmodel       <- fit0@Model
  lavsamplestats <- fit0@SampleStats
  lavdata        <- fit0@Data
  lavoptions     <- fit0@Options
  lavpartable    <- fit0@ParTable
  n              <- fit0@SampleStats@ntotal

  # Starting values
  start <- lavaan::lav_model_get_parameters(lavmodel)  # this will be packed

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
      x = start,
      lower = lower,
      upper = upper,
      n = n,

      # Information about constraints
      lavmodel_eq.constraints = lavmodel@eq.constraints,
      lavmodel_eq.constraints.k0 = lavmodel@eq.constraints.k0,
      lavmodel_eq.constraints.K = lavmodel@eq.constraints.K,

      lavmodel_ceq.simple.only = lavmodel@ceq.simple.only,
      lavmodel_ceq.simple.K = lavmodel@ceq.simple.K,

      # Information about the model
      loglik = loglik(start, lavmodel, lavsamplestats, lavdata, lavoptions,
                      bias_reduction = isFALSE(is_ML | is_eRBM),
                      plugin_pen = plugin_pen, bounds = list(lower, upper),
                      verbose = FALSE, fn.scale = fn.scale),


      grad_loglik = grad_loglik(start, lavmodel, lavsamplestats, lavdata,
                                lavoptions),
      j = information_matrix(start, lavmodel, lavsamplestats, lavdata,
                             lavoptions, kind = information),
      jinv = try(solve(
        information_matrix(start, lavmodel, lavsamplestats, lavdata,
                           lavoptions, kind = information),
      )),
      e = information_matrix(start, lavmodel, lavsamplestats, lavdata,
                             lavoptions, kind = "first.order"),
      penalty = penalty(start, lavmodel, lavsamplestats, lavdata, lavoptions),
      bias = bias(start, lavmodel, lavsamplestats, lavdata, lavoptions),
      information = information
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
      x = x,
      ...,
      bias_reduction = isFALSE(is_ML | is_eRBM),
      plugin_pen = plugin_pen,
      bounds = list(lower = lower, upper = upper),
      verbose = verbose,
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
    b <- bias(
      x = res$par,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
  est <- res$par - b

  elapsed_time <- proc.time() - start_time
  elapsed_time <- elapsed_time["elapsed"]

  # Implied model covariance matrix --------------------------------------------
  lavmodel <- lavaan::lav_model_set_parameters(lavmodel, est)
  Sigma    <- lavaan::lav_model_implied(lavmodel)$cov[[1]]

  # Get max grad ---------------------------------------------------------------
  if (isTRUE(maxgrad)) {
    hinv <- try(solve(numDeriv::hessian(
      func = obj_fun,
      x = res$par,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )))
    max_scores <- numDeriv::grad(
      func = obj_fun,
      x = res$par,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
    if (inherits(hinv, "try-error")) {
      scaled_grad <- NA
    } else {
      scaled_grad <- as.numeric(hinv %*% max_scores)
    }
  } else {
    scaled_grad <- NULL
  }

  # Standard errors ------------------------------------------------------------
  j <- information_matrix(
    x = as.numeric(est),
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    kind = information
  )
  if (lavmodel@ceq.simple.only) {
    K <- lavmodel@ceq.simple.K
    j <- t(K) %*% j %*% K
  }
  if (check_mat(j)) {
    sds <- rep(NA, length(est))
    jinv <- NULL
  } else {
    jinv <- try(solve(j), silent = !TRUE)
    sds <- sqrt(diag(jinv))
  }

  # Unpack estimators and se
  idx <- lavpartable$free[lavpartable$free > 0]
  est <- est[idx]
  sds <- sds[idx]
  names(est) <- names(coef(fit0))

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
