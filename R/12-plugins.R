pen_ridge <- function(x, call = FALSE, ...) {
  if (isTRUE(call)) return("Ridge penalty")
  sum(x ^ 2)
}

pen_ridge_bound <- function(x, lb, ub, call = FALSE, ...) {

  if (isTRUE(call)) return("Bounded ridge penalty")

  fb <- which(!(is.infinite(lb) | is.infinite(ub)))

  below_penalty <- ifelse(x < lb, (lb - x) ^ 2, 0)
  above_penalty <- ifelse(x > ub, (x - ub) ^ 2, 0)

  sum(below_penalty[fb] + above_penalty[fb])
}

# From IK 7/12/24
#
# huber <- function(x, thres = 1) {
#   ifelse(abs(x) < thres, -0.5 * x^2, thres * (thres / 2 - abs(x)))
# }
#
# huber_doubly_bounded <- function(x, low, upp, thres = 1) {
#   x <- log((x - low) / (upp - x))
#   huber(x, thres = thres)
# }
#
# huber_positive <- function(x, low = 0, thres = 1) {
#   x <- log(x - low)
#   huber(x, thres = thres)
# }

pen_huber <- function(x, lb, ub, thres = 1, call = FALSE, ...) {

  if (isTRUE(call)) return("Huber penalty")

  # A bit arbitrary but force values of x outside bounds to be equal to bounds
  x[x > ub] <- ub[x > ub]
  x[x < lb] <- lb[x < lb]

  y <-
    ifelse(
      is.finite(lb),
      ifelse(
        is.finite(ub),
        log((x - lb) / (ub - x)),  # both lb and ub finite
        log(x - lb)  # only lb finite
      ),
      0  # both lb and ub infinite, no penalty
    )
  # return(data.frame(x = x, lb = lb, ub = ub, y = y))
  # y <- y[is.finite(y)]
  sum(ifelse(abs(y) < thres, 0.5 * y ^ 2, -thres * (thres / 2 - abs(y))))
}

# x <- c(275, 0, 50, 0, 20, 1000)
# lb <- c(0, -Inf, 0, -Inf, -10694, -40)
# ub <- c(10694, Inf, 10694, Inf, 10694, 953)
# pen_huber(x, lb, ub)
