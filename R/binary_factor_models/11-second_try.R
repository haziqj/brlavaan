library(tidyverse)
library(lavaan)
library(lavaan.bingof)
# pak::pak("haziqj/lavaan.bingof")
library(numDeriv)

model_no <- 1
mod <- txt_mod(model_no)
dat <- gen_data_bin(model_no)
theta0 <- get_true_values(model_no)

# Using lavaan
fit_lav <- sem(mod, dat, std.lv = TRUE, estimator = "PML")
# theta_hat <- parTable(fit_lav)$est[1:10]
theta_hat <- coef(fit_lav)

# Write the J and H functions
JJJ <- function(theta) {
  # Fix param values
  # fit_lav_fxd <- sem(mod, dat, std.lv = TRUE, estimator = "PML",
  #                    start = parTable(fit_lav)$est[1:10], baseline = FALSE,
  #                    optim.method = "none")
  fit_lav_fxd <- sem(mod, dat, std.lv = TRUE, estimator = "PML",
                     do.fit = FALSE)
  # extract slots
  lavmodel       <- fit_lav_fxd@Model
  lavsamplestats <- fit_lav_fxd@SampleStats
  lavdata        <- fit_lav_fxd@Data
  lavoptions     <- fit_lav_fxd@Options

  # update lavmodel with 'new' set of values
  lavmodel <- lav_model_set_parameters(lavmodel, x = theta)
  lavaan:::lav_model_information_firstorder(
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
}

HHH <- function(theta) {
  # Fix param values
  # fit_lav_fxd <- sem(mod, dat, std.lv = TRUE, estimator = "PML",
  #                    start = parTable(fit_lav)$est[1:10], baseline = FALSE,
  #                    optim.method = "none")
  fit_lav_fxd <- sem(mod, dat, std.lv = TRUE, estimator = "PML",
                     do.fit = FALSE)
  # extract slots
  lavmodel       <- fit_lav_fxd@Model
  lavsamplestats <- fit_lav_fxd@SampleStats
  lavdata        <- fit_lav_fxd@Data
  lavoptions     <- fit_lav_fxd@Options
  lavcache       <- fit_lav_fxd@Cache
  lavimplied     <- fit_lav_fxd@implied
  lavh1          <- fit_lav_fxd@h1

  # update lavmodel with 'new' set of values
  lavmodel <- lav_model_set_parameters(lavmodel, x = theta)

  # lavaan:::lav_model_information(
  #   lavmodel       = lavmodel,
  #   lavsamplestats = lavsamplestats,
  #   lavdata        = lavdata,
  #   lavcache       = lavcache,
  #   lavimplied     = lavimplied,
  #   lavh1          = lavh1,
  #   lavoptions     = lavoptions,
  #   extra          = FALSE,
  #   augmented      = TRUE,
  #   inverted       = TRUE,
  #   use.ginv       = FALSE)
  H <- lavaan:::lav_model_hessian(
    lavmodel       = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata        = lavdata,
    lavoptions     = lavoptions,
    lavcache       = lavcache,
    ceq.simple     = FALSE
  )
  H / lavsamplestats@ntotal  # UNIT INFORMATION
}

AAA <- function(.theta) {
  tmp <- function(x) {
    Hinv <- HHH(x) |> MASS::ginv()
    J    <- JJJ(x)
    -0.5 * sum(diag(Hinv %*% J))
  }
  numDeriv::grad(tmp, .theta)
}

# Bias reduction
(A <- AAA(theta_hat))
Hinv <- lavaan.bingof:::get_sensitivity_inv_mat(fit_lav)
theta_tilde <- theta_hat + Hinv %*% A
theta_tilde
