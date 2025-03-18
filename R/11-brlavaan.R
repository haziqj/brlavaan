#' Fit a latent variable model using bias-reducing methods
#'
#' Fit a latent variable model using bias-reducing methods
#'
#' The `pen_ridge()` function applies a ridge regression style penalty $f(x) =
#' || x ||^2$ that shrinks the parameters to zero. The `pen_ridge_bound()`
#' function applies a penalty that shrinks the parameters to the bounds of the
#' parameter space. The bounds are calculated by `{lavaan}` -- see the paper for
#' more details.
#'
#' The `info_pen` and `info_bias` arguments default to `"observed"`. It is not
#' recommended to change this to say `"expected"`, especially for small sample
#' sizes as the bias-reducing properties of the estimators are not guaranteed.
#'
#' @param model A description of the user-specified model. Typically, the model
#'   is described using the lavaan model syntax. See [model.syntax] for more
#'   information. Alternatively, a parameter table (eg. the output of the
#'   [lavaanify()] function) is also accepted.
#' @param data An optional data frame containing the observed variables used in
#'   the model. If some variables are declared as ordered factors, lavaan will
#'   treat them as ordinal variables.
#' @param estimator The estimator to use. **Currently only "ML" is supported.**
#' @param estimator.args A list containing RBM arguments. Possible arguments are
#'   - `rbm`: The type of RBM method to use. One of `"none"`, `"explicit"`, or
#'   `"implicit"` (although, short forms are accepted, e.g. `"exp"`, `"iRBM"`,
#'   etc.)
#'   - `plugin_pen`: The type of penalty to use. One of `NULL`, `"pen_ridge"`,
#'   or `"pen_ridge_bound"`.
#'   - `info_pen` The type of information matrix to use for the penalty term.
#'   - `info_bias` The type of information matrix to use for the bias term of
#'   the explicit reduced bias method.
#' @param information The type of information matrix to use for calculation of
#'   standard errors. Defaults to `"observed"`, although `"expected"` and
#'   `"first.order"` is also permitted.
#' @param lavfun The lavaan function to use. Default is "sem".
#' @param ... Additional arguments to pass to the [lavaan] function.
#'
#' @return An object of class `brlavaan` which is a subclass of the
#'   [lavaan-class] class.
#' @export
brlavaan <- function(
    model,
    data,
    estimator = "ML",
    estimator.args = list(rbm = "implicit", plugin_pen = NULL),
    information = "observed",
    lavfun = "sem",
    ...
) {

  if (length(estimator.args$info_pen) > 0) {
    if (estimator.args$info_pen == "expected")
      cli::cli_alert_warning("Using expected information matrix for the penalty term is not guaranteed to produce bias-reducing properties in small samples.")
  }

  fit <- fit_sem(
    model = model,
    data = data,
    debug = FALSE,
    lavfun = lavfun,
    estimator = estimator,
    estimator.args = estimator.args,
    information = information,
    ...
  )

  out <- create_lav_from_fitsem(fit, model, data, ...)
  new("brlavaan", out)
}

#' Fit a structural equation model using empirical bias-reducing methods
#'
#' @inherit brlavaan params return
#' @export
brsem <- function(
    model,
    data,
    estimator = "ML",
    estimator.args = list(rbm = "implicit", plugin_pen = NULL),
    information = "observed",
    ...
) {

  brlavaan(
    model = model,
    data = data,
    estimator = estimator,
    estimator.args = estimator.args,
    information = information,
    lavfun = "sem",
    ...
  )
}

#' Fit confirmatory factor analysis model using empirical bias-reducing methods
#'
#' @inherit brlavaan params return
#' @export
brcfa <- function(
    model,
    data,
    estimator = "ML",
    estimator.args = list(rbm = "implicit", plugin_pen = NULL),
    information = "observed",
    ...
) {

  brlavaan(
    model = model,
    data = data,
    estimator = estimator,
    estimator.args = estimator.args,
    information = information,
    lavfun = "cfa",
    ...
  )
}

#' Fit growth curve models using empirical bias-reducing methods
#'
#' @inherit brlavaan params return
#' @export
brgrowth <- function(
    model,
    data,
    estimator = "ML",
    estimator.args = list(rbm = "implicit", plugin_pen = NULL),
    information = "observed",
    ...
  ) {

  brlavaan(
    model = model,
    data = data,
    estimator = estimator,
    estimator.args = estimator.args,
    information = information,
    lavfun = "growth",
    ...
  )
}

create_lav_from_fitsem <- function(
    fit,
    model,
    data,
    ...
  ) {

  # Prepare output (blank lavaan object)
  x <- fit$coefficients
  lavargs <- list(...)
  lavargs$model <- model
  lavargs$data <- data
  lavargs$do.fit <- FALSE
  lavargs$method <- NULL  # if using old ways of specifying RBM method
  lavargs$information <- fit$information
  lavargs$start <- x
  lavargs$estimator.args <- list(
    rbm = fit$rbm,
    plugin_pen = fit$plugin_pen
  )
  fit0 <- do.call(get(fit$lavfun, envir = asNamespace("lavaan")), lavargs)

  # Change version slot
  fit0@version <- as.character(packageVersion("brlavaan"))

  # Change timing slot
  fit0@timing$optim <- fit0@timing$optim + fit$timing
  fit0@timing$total <- fit0@timing$total + fit$timing

  # Change Model and implied slots
  fit0@Model <- lavaan::lav_model_set_parameters(fit0@Model, x)
  fit0@implied <- lavaan::lav_model_implied(fit0@Model)

  # Change ParTable and pta slots
  pt <- lavaan::partable(fit0)
  pt$est[pt$free > 0] <- x
  pt$se <- 0
  pt$se[pt$free > 0] <- fit$stderr
  fit0@ParTable <- as.list(pt)
  fit0@pta$names <- names(pt)

  # Change Options slot
  fit0@Options$estimator <- fit$estimator
  # fit0@Options$estimator.args <- list(method = "eRBM")
  # fit0@Options$test <- "standard"
  fit0@Options$se <- "standard"
  fit0@Options$do.fit <- TRUE

  # Change Fit slot
  fit0@Fit@x <- x
  fit0@Fit@se <- fit$stderr
  fit0@Fit@iterations <- fit$optim$iterations
  fit0@Fit@converged <- fit$optim$convergence == 0L

  # Change optim slot
  fit0@optim$x <- x
  # fit0@optim$dx <- 0
  fit0@optim$npar <- length(x)
  fit0@optim$fx <- fit0@Fit@fx
  fit0@optim$fx.group <- fit0@Fit@fx.group
  fit0@optim$iterations <- fit$optim$iterations
  fit0@optim$converged <- fit$optim$convergence == 0L

  # Change loglik slot
  # fit0@loglik$estimator <-
  #   if (fit$estimator == "ML") "ML"
  #   else if (fit$estimator == "IBRM") "IMP-BR ML"
  #   else if (fit$estimator == "IBRMP") "IMP-BR ML"
  #   else if (fit$estimator == "EBRM") "EXP-BR ML"
  # Change vcov slot
  fit0@vcov$se <- "standard"
  fit0@vcov$vcov <- fit$vcov

  # fit0@test <- fit_lav@test
  # fit0@baseline <- fit_lav@baseline

  # Include the entire output of fit_sem
  fit0@external <- fit

  fit0
}

