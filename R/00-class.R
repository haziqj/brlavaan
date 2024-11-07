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

setMethod("show", "brlavaan", function(object) {
  cat("br")
  callNextMethod()
})

setMethod("summary", "brlavaan", function(object, ...) {
  class(object) <- "lavaan"
  cat("br")
  out <- capture.output(callNextMethod())

  if (object@Options$estimator == "IBRMP")
    penalty_report <-
      "  Penalty                                         TRUE"
  else
    penalty_report <-
      "  Penalty                                        FALSE"

  out <- append(out, penalty_report, after = 3)
  cat(out, sep = "\n")

})

setMethod("coef", "brlavaan", function(object, ...) {
  class(object) <- "lavaan"
  callNextMethod()
})
