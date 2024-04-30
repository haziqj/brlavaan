library(lavaan)
library(tidyverse)
library(furrr)
plan(multisession, workers = parallel::detectCores() - 2)
source("R/factor_models/10-extract_H_and_J_matrices.R")
source("R/factor_models/20-jackknife_bootstrap.R")
B <- 1000  # no of sims

growth_sim_fn <- function(i = 1, n = 200) {
  # Generate data
  mod <- "
  i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4 + 1*y5 + 1*y6 + 1*y7 + 1*y8 + 1*y9 + 1*y10
  s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4 + 4*y5 + 5*y6 + 6*y7 + 7*y8 + 8*y9 + 9*y10

  i ~~ 550 * i
  s ~~ 100 * s
  i ~~ 40 * s
  y1 ~~ 500 * y1
  y2 ~~ 500 * y2
  y3 ~~ 500 * y3
  y4 ~~ 500 * y4
  y5 ~~ 500 * y5
  y6 ~~ 500 * y6
  y7 ~~ 500 * y7
  y8 ~~ 500 * y8
  y9 ~~ 500 * y9
  y10 ~~ 500 * y10
  "
  dat <- simulateData(model = mod, sample.nobs = n)

  fit <- sem(
    model = "
    i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4 + 1*y5 + 1*y6 + 1*y7 + 1*y8 + 1*y9 + 1*y10
    s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4 + 4*y5 + 5*y6 + 6*y7 + 7*y8 + 8*y9 + 9*y10
    ",
    data = dat
  )

  tibble(
    i     = i,
    par   = names(coef(fit)),
    truth = c(0.8, 0.6, 1, 1, 1, 1),  # CHANGE
    ml    = coef(fit),
    ebrm  = rb_empr(fit),
    jack  = rb_jack(fit, dat),
    boot  = rb_boot(fit, dat)
  )
}

res <-
  furrr::future_map(
    .x = 1:B,
    .f = possibly(sim_fn, NA),
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
  )
