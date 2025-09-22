# Toy example for N(0, sigma2) to show effect of rbm
library(tidyverse)

set.seed(123)
n      <- 50
sigma2 <- 1  # theta
x      <- rnorm(n, mean = 0, sd = sqrt(sigma2))
S      <- sum(x^2)

# Ingredients
loglik <- function(theta, xx = x) sum(dnorm(xx, sd = sqrt(theta), log = TRUE))
score <- function(theta, xx = x) {
  sum(-1 / (2 * theta) + (xx - 0) ^ 2 / (2 * theta ^ 2))
}


e <- function(theta) {
  tmp <- -1 / (2 * theta) + (x - 0) ^ 2 / (2 * theta ^ 2)
  sum(tmp ^ 1)
}

j <- function(theta) {
  tmp <- 1 / (2 * theta ^ 2) - (x - 0) ^ 2 / (theta ^ 3)
  -1 * sum(tmp)
}

# enum <- function(theta) {
#   res <- x
#   for (i in seq_along(res)) {
#     res[i] <- numDeriv::grad(function(theta) dnorm(x[i], sd = sqrt(theta), log = TRUE), theta)
#   }
#   sum(res ^ 2)
# }
#
# jnum <- function(theta) {
#   as.numeric(-1 * numDeriv::hessian(loglik, theta))
# }
#
#
penalty <- function(theta) -0.5 * e(theta) / j(theta)
#
# theta_hat <- nlminb(1, function(x) -1 * loglik(x), lower = 0.01, upper = 5)$par
# theta_til <- nlminb(1, function(x) -(loglik(x) + penalty(x)), lower = 0.01, upper = 5)$par
# c(theta_hat, theta_til, sigma2)
#
# tibble(theta = seq(0.3, 1.5, length = 100)) |>
#   mutate(
#     # loglik = map_dbl(theta, loglik),
#     pen = map_dbl(theta, penalty),
#     # pen = loglik + pen
#   ) |>
#   # pivot_longer(cols = c(loglik, pen), names_to = "which", values_to = "value") |>
#   ggplot(aes(theta, pen)) +
#   geom_line() +
#   theme_bw()
#
# tibble(theta = seq(0.5, 1.5, length = 100)) |>
#   rowwise() |>
#   mutate(
#     e = e(theta),
#     j = j(theta),
#     jinv = 1 / j,
#     pen = penalty(theta)
#   ) |>
#   pivot_longer(cols = c(e), names_to = "which", values_to = "value") |>
#   ggplot(aes(theta, value, col = which)) +
#   geom_line() +
#   theme_bw()


## -----------------------------------------------------------------------------
set.seed(123)
n <- 15
B <- 50
sigma2 <- 1  # theta
X <- lapply(1:B, function(i) {
  rnorm(n, mean = 0, sd = sqrt(sigma2))
})

bump <- 3 / n # exagerate the bias, diminishes as n -> Inf

plot_df <-
  tibble(b = 1:B) |>
  mutate(
    X = map(b, \(i) X[[i]]),
    theta = list(seq(0.32, 1.25, length = 1000)),
    score = map2(X, theta, \(xx, tt) 2 * map_dbl(tt + bump, score, xx))
  ) |>
  select(b, theta, score) |>
  unnest(c(theta, score))

plot_df2 <-  summarise(plot_df, score = mean(score), .by = theta)
theta_hat <-
  plot_df2 |>
  mutate(absscore = abs(score)) |>
  filter(absscore == min(absscore)) |>
  slice(1) |>
  pull(theta)
bias <- theta_hat - sigma2
plot_df3 <-
  plot_df2 |>
  mutate(theta = theta - bias) |>
  filter(theta <= 1.25)

ggplot(plot_df, aes(theta, score, group = b)) +
  geom_line(linewidth = 0.1, col = "gray40") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(
    data = plot_df2,
    aes(theta, score),
    col = "blue",
    inherit.aes = FALSE,
    linewidth = 1
  ) +
  geom_line(
    data = plot_df3,
    aes(theta, score),
    inherit.aes = FALSE,
    col = "red3",
    linewidth = 1
  ) +
  scale_x_continuous(
    breaks = c(theta_hat, sigma2),
    labels = c(expression(hat(theta), theta[0])),
    name = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  coord_cartesian(ylim = c(-8, 28)) +
  labs(y = expression(Score~U(theta)))


