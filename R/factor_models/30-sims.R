library(lavaan)
library(tidyverse)
library(furrr)
plan(multisession, workers = parallel::detectCores() - 2)
source("R/factor_models/10-extract_H_and_J_matrices.R")
source("R/factor_models/20-jackknife_bootstrap.R")
B <- 1000  # no of sims

sim_fn <- function(i = 1) {
  # Generate data
  mod <- "fx =~ 1*x1 + 0.8*x2 + 0.6*x3"
  n   <- 200
  dat <- simulateData(model = mod, sample.nobs = n)
  fit <- cfa(model = "fx =~ x1 + x2 + x3", data = dat)

  tibble(
    i     = i,
    par   = names(coef(fit)),
    truth = c(0.8, 0.6, 1, 1, 1, 1),
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

save(res, file = "R/factor_models/sims_toyexample.RData")

# Mean bias
do.call(rbind, res) |>
  summarise(
    across(c(ml, ebrm, jack, boot), \(x) mean(x - truth, na.rm = TRUE)),
    .by = par
  )
