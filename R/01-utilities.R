validate_rbm <- function(x) {
  valid_rbm <- c("none", "explicit", "implicit")

  if (length(x) == 0) x <- "none"
  x[is.na(x)] <- "none"
  x[isFALSE(as.logical(x))] <- "none"

  x[startsWith(x, "e")] <- "explicit"
  x[startsWith(x, "i")] <- "implicit"

  rlang::arg_match(x, valid_rbm)
}

get_lav_stuff <- function(fit) {
  # Utility function to extract lavaan stuff
  list(
    lavmodel       = fit@Model,
    lavsamplestats = fit@SampleStats,
    lavdata        = fit@Data,
    lavoptions     = fit@Options
  )
}

get_estimator.args <- function(x) {
  if (inherits(x, "brlavaan")) {
    rbm <- x@Options$estimator.args$rbm
    plugin_pen <- x@Options$estimator.args$plugin_pen
  } else {
    rbm <- x$rbm
    plugin_pen <- x$plugin_pen
  }
  list(rbm = rbm, plugin_pen = plugin_pen)
}

# To remove R CMD CHECK NOTE "no visible binding"
globalVariables(c("rbm", "plugin_pen"))

#' Predicates for `brlavaan` objects
#'
#' @param x An object of class `brlavaan`, or an output from [`fit_sem()`].
#' @param quietly If `TRUE`, suppresses messages.
#'
#' @return `TRUE` or `FALSE`, or if not suppressed, a message indicating which
#'   plugin penalty was used.
#' @name predicates
NULL

#' @rdname predicates
#' @export
is_ML <- function(x) {
  list2env(get_estimator.args(x), environment())
  rbm == "none"
}

#' @rdname predicates
#' @export
is_eRBM <- function(x) {
  list2env(get_estimator.args(x), environment())
  rbm == "explicit"
}

#' @rdname predicates
#' @export
is_iRBM <- function(x) {
  list2env(get_estimator.args(x), environment())
  rbm == "implicit" & is.null(plugin_pen)
}

#' @rdname predicates
#' @export
is_iRBMp <- function(x, quietly = FALSE) {
  list2env(get_estimator.args(x), environment())

  out <- rbm == "implicit" & !is.null(plugin_pen)

  rbm <- paste0(
    toupper(substr(rbm, 1, 1)),
    tolower(substr(rbm, 2, nchar(rbm)))
  )
  plugin_pen <- plugin_pen(call = TRUE)

  if (isFALSE(quietly))
    cli::cli_alert_info("is_iRBMp: RBM = {rbm}, Plugin penalty = {plugin_pen}")

  attr(out, "plugin_pen") <- plugin_pen
  out
}

# # I thought it would be useful to have functions for each kind of information
#
# # Observed information
# information_observed <- function(
#     theta,
#     lavmodel,
#     lavsamplestats,
#     lavdata,
#     lavoptions
#   ) {
#   information_matrix(
#     theta = theta,
#     lavmodel = lavmodel,
#     lavsamplestats = lavsamplestats,
#     lavdata = lavdata,
#     lavoptions = lavoptions,
#     kind = "observed"
#   )
# }
#
# # Expected information
# information_expected <- function(
    #     theta,
#     lavmodel,
#     lavsamplestats,
#     lavdata,
#     lavoptions
#   ) {
#   information_matrix(
#     theta = theta,
#     lavmodel = lavmodel,
#     lavsamplestats = lavsamplestats,
#     lavdata = lavdata,
#     lavoptions = lavoptions,
#     kind = "expected"
#   )
# }
#
# # Outer product of casewise scores
# information_firstorder <- function(
    #     theta,
#     lavmodel,
#     lavsamplestats,
#     lavdata,
#     lavoptions
#   ) {
#   information_matrix(
#     theta = theta,
#     lavmodel = lavmodel,
#     lavsamplestats = lavsamplestats,
#     lavdata = lavdata,
#     lavoptions = lavoptions,
#     kind = "first.order"
#   )
# }
