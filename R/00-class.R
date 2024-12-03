#' brlavaan Class
#'
#' This is a class that extends the lavaan class.
#'
#' @importFrom lavaan lavaan
#' @export
setClass(
  Class = "brlavaan",
  contains = "lavaan"
)

rbmpen_report <- function(object) {
  if (is_ML(object)) {
    rbm_report <-
      "  Bias reduction method                           NONE"
    penalty_report <-
      "  Plugin penalty                                  NONE"
  } else if (is_eRBM(object)) {
    rbm_report <-
      "  Bias reduction method                       EXPLICIT"
    penalty_report <-
      "  Plugin penalty                                  NONE"
  } else if (is_iRBM(object)) {
    rbm_report <-
      "  Bias reduction method                       IMPLICIT"
    penalty_report <-
      "  Plugin penalty                                  NONE"
  } else {
    tmp <- is_iRBMp(object, quietly = TRUE)
    plugin_penalty <- attr(tmp, "plugin_penalty")
    rbm_report <-
      "  Bias reduction method                       IMPLICIT"
    penalty_report <-
      str_pad(plugin_penalty, 54 - nchar("  Plugin penalty "), side = "left")
    penalty_report <- sprintf("  Plugin penalty %s", penalty_report)
  }

  list(rbm_report = rbm_report, penalty_report = penalty_report)
}

setMethod("show", "brlavaan", function(object) {

  rbm_report <- rbmpen_report(object)$rbm_report
  penalty_report <- rbmpen_report(object)$penalty_report

  class(object) <- "lavaan"
  cat("br")
  out <- capture.output(callNextMethod())
  out <- out[seq_len(length(out) / 2)]  # not sure why there are duplicates
  out <- append(out, rbm_report, after = 3)
  out <- append(out, penalty_report, after = 4)
  cat(out, sep = "\n")

})

setMethod("summary", "brlavaan", function(object, ...) {

  rbm_report <- rbmpen_report(object)$rbm_report
  penalty_report <- rbmpen_report(object)$penalty_report

  class(object) <- "lavaan"
  cat("br")
  out <- capture.output(callNextMethod())
  out <- append(out, rbm_report, after = 3)
  out <- append(out, penalty_report, after = 4)
  cat(out, sep = "\n")

})

setMethod("coef", "brlavaan", function(object, ...) {
  class(object) <- "lavaan"
  callNextMethod()
})
