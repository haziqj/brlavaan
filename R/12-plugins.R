pen_ridge <- function(x, lb, ub) {
  sum(x ^ 2)
}

pen_ridge_bound <- function(x, lb, ub) {
  fb <- which(!(is.infinite(lb) | is.infinite(ub)))

  below_penalty <- ifelse(x < lb, (lb - x) ^ 2, 0)
  above_penalty <- ifelse(x > ub, (x - ub) ^ 2, 0)

  sum(below_penalty[fb] + above_penalty[fb])
}
