library(tidyverse)
library(brlavaan)

## ----- Generate data set -----------------------------------------------------
set.seed(1234)
X <- rnorm(100)
M <- 0.5 * X + rnorm(100)
Y <- 0.7 * M + rnorm(100)

dat <- tibble(X = X, Y = Y, M = M)
model <- '
  # direct effect
  Y ~ c*X
  # mediator
  M ~ a*X
  Y ~ b*M
  # indirect effect (a*b)
  ab := a*b
  # total effect
  total := c + (a*b)
'

## ----- Fit models ------------------------------------------------------------
fit <- brsem(model, data = dat)
summary(fit)
