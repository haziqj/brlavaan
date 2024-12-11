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
    plugin_penalty <- x@Options$estimator.args$plugin_penalty
  } else {
    rbm <- x$rbm
    plugin_penalty <- x$plugin_penalty
  }
  list(rbm = rbm, plugin_penalty = plugin_penalty)
}

# To remove R CMD CHECK NOTE "no visible binding"
globalVariables(c("rbm", "plugin_penalty"))

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
  isFALSE(rbm)
}

#' @rdname predicates
#' @export
is_eRBM <- function(x) {
  list2env(get_estimator.args(x), environment())
  rbm == "eRBM"
}

#' @rdname predicates
#' @export
is_iRBM <- function(x) {
  list2env(get_estimator.args(x), environment())
  rbm == "iRBM" & is.null(plugin_penalty)
}

#' @rdname predicates
#' @export
is_iRBMp <- function(x, quietly = FALSE) {
  list2env(get_estimator.args(x), environment())
  plugin_penalty <- plugin_penalty(call = TRUE)

  if (isFALSE(quietly))
    cli::cli_alert_info("is_iRBMp: rbm = {rbm}, plugin_penalty = {plugin_penalty}")
  out <- rbm == "iRBM" & !is.null(plugin_penalty)
  attr(out, "plugin_penalty") <- plugin_penalty
  out
}
