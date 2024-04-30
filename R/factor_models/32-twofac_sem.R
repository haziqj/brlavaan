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
  eta1 =~ 1*y1 + 0.7*y2 + 0.6*y3
  eta2 =~ 1*y4 + 0.7*y5 + 0.6*y6
  eta2 ~ 0.25*eta1

  eta1 ~~ 1*eta1
  eta2 ~~ 1*eta2
  y1 ~~ 0.25*y1
  y4 ~~ 0.25*y4
  y2 ~~ 0.09*y2
  y5 ~~ 0.09*y5
  y3 ~~ 0.1225*y3
  y6 ~~ 0.1225*y6
  "
  dat <- simulateData(model = mod, sample.nobs = n)

  fit <- sem(
    model = "
    eta1 =~ y1 + y2 + y3
    eta2 =~ y4 + y5 + y6
    eta2 ~ eta1
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
