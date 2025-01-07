# A latent interaction model examines the interaction between two latent
# variables. For example, we might study how Motivation (M) moderates the effect
# of Skill (S) on Performance (P) in a workplace or educational setting.
library(tidyverse)
library(brlavaan)
library(semTools)

## ----- Generate data set -----------------------------------------------------
n <- 5000  # sample size

# Define factor loadings
lambda_sk <- c(0.8, 0.7, 0.6)  # skill loadings
lambda_mo <- c(0.9, 0.8, 0.7)  # motivation loadings
lambda_pe <- c(0.85, 0.75, 0.65)  # performance loadings

# Correlation between sk and mo
cor_sk_mo <- 0.5

# Generate latent variables
eta <- mvtnorm::rmvnorm(n, sigma =  matrix(c(1, 0.5, 0.5, 1), nrow = 2))
sk <- eta[, 1]
mo <- eta[, 2]

# Define regression coefficients for pe ~ sk + mo + sk*mo
beta_sk <- 1.2
beta_mo <- 0.8
beta_skmo <- 0.3
pe <- beta_sk * sk + beta_mo * mo + beta_skmo * sk * mo + rnorm(n, sd = 0.5)

# Generate indicators for sk
sk1 <- lambda_sk[1] * sk + rnorm(n, sd = sqrt(1 - lambda_sk[1] ^ 2))
sk2 <- lambda_sk[2] * sk + rnorm(n, sd = sqrt(1 - lambda_sk[2] ^ 2))
sk3 <- lambda_sk[3] * sk + rnorm(n, sd = sqrt(1 - lambda_sk[3] ^ 2))

# Generate indicators for mo
mo1 <- lambda_mo[1] * mo + rnorm(n, sd = sqrt(1 - lambda_mo[1] ^ 2))
mo2 <- lambda_mo[2] * mo + rnorm(n, sd = sqrt(1 - lambda_mo[2] ^ 2))
mo3 <- lambda_mo[3] * mo + rnorm(n, sd = sqrt(1 - lambda_mo[3] ^ 2))

# Generate indicators for pe
pe1 <- lambda_pe[1] * pe + rnorm(n, sd = sqrt(1 - lambda_pe[1] ^ 2))
pe2 <- lambda_pe[2] * pe + rnorm(n, sd = sqrt(1 - lambda_pe[2] ^ 2))
pe3 <- lambda_pe[3] * pe + rnorm(n, sd = sqrt(1 - lambda_pe[3] ^ 2))

# Combine data into a data frame
dat <- tibble(sk1, sk2, sk3, mo1, mo2, mo3, pe1, pe2, pe3)

## ----- Fit models ------------------------------------------------------------

# Model fit without interaction
mod <- "
  ski =~ sk1 + sk2 + sk3
  mot =~ mo1 + mo2 + mo3
  per =~ pe1 + pe2 + pe3

  per ~ ski + mot
"
fit0 <- sem(mod, dat)

# Model fit with interaction
dat_dmc <- as_tibble(indProd(
  data = dat,
  var1 = c("sk1", "sk2", "sk3"),
  var2 = c("mo1", "mo2", "mo3"),
  match = FALSE,  # all pairs strategy
  meanC = FALSE,  # only one of meanC, residualC, doubleMC should be TRUE
  residualC = FALSE,
  doubleMC = TRUE
))

mod_int <- "
  # Measurement model
  ski =~ sk1 + sk2 + sk3
  mot =~ mo1 + mo2 + mo3
  per =~ pe1 + pe2 + pe3
  int =~ sk1.mo1 + sk1.mo2 + sk1.mo3 +
         sk2.mo1 + sk2.mo2 + sk2.mo3 +
         sk3.mo1 + sk3.mo2 + sk3.mo3

  # Structural regression
  per ~ ski + mot + int

  # Residual covariances with equality constraints between same items
  sk1.mo1 ~~ v1 * sk1.mo2 + v1 * sk1.mo3
  sk1.mo2 ~~ v1 * sk1.mo3
  sk2.mo1 ~~ v2 * sk2.mo2 + v2 * sk2.mo3
  sk2.mo2 ~~ v2 * sk3.mo3
  sk3.mo1 ~~ v3 * sk3.mo2 + v3 * sk3.mo3
  sk3.mo2 ~~ v3 * sk3.mo3
  sk1.mo1 ~~ v4 * sk2.mo1 + v4 * sk3.mo1
  sk2.mo1 ~~ v4 * sk3.mo1
  sk1.mo2 ~~ v5 * sk2.mo2 + v5 * sk3.mo3
  sk2.mo2 ~~ v5 * sk3.mo2
  sk1.mo3 ~~ v6 * sk2.mo3 + v6 * sk3.mo3
  sk2.mo3 ~~ v6 * sk3.mo3
"

fit1 <- sem(mod_int, data = dat_dmc, std.lv = TRUE)
fit2 <- brsem(mod_int, data = dat_dmc, method = "eRBM", std.lv = TRUE, meanstructure = TRUE)
fit3 <- brsem(mod_int, data = dat_dmc, method = "iRBM", std.lv = TRUE, meanstructure = TRUE)

# Simple slopes
probe <- probe2WayMC(fit3, c("ski", "mot", "int"), "per", "mot", c(-1, 0, 1))
plotProbe(probe, c(-3, 3))

## ----- Compare models --------------------------------------------------------
imap(
  .x = list(ML = fit1, eRBM = fit2, iRBM = fit3),
  .f = \(x, idx) {
  partable(x) |>
    as_tibble() |>
    dplyr::filter(op == "~") |>
    dplyr::select(lhs, op, rhs, est, se) |>
    unite(col = "coef", lhs, op, rhs, sep = "") |>
    mutate(
      truth = c(1.2, 0.8, 0.3),
      bias = est - truth,
      method = idx
    )
}) |>
  do.call(what = rbind) |>
  pivot_wider(id_cols = coef, names_from = method, values_from = c(bias))

# Some thoughts: I'm not sure the comparison is valid. Because the true data
# generating mechanism comes from a "true" latent interaction model. But the
# estimation is based on a "product-indicator" approach. So the parameters will
# not be correctly estimated.



