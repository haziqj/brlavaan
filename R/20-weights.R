library(tidyverse)
library(lavaan)
library(lavaan.bingof)

model_no <- 1
mod <- txt_mod(model_no)
dat <- gen_data_bin_wt(model_no, n = 1000)
theta0 <- get_true_values(model_no)

# Using lavaan
fit_lav <- sem(mod, dat, std.lv = TRUE, estimator = "PML")
# theta_hat <- parTable(fit_lav)$est[1:10]
theta_hat <- coef(fit_lav)

# Bias reduction
(A <- AAA(theta_hat))
Hinv <- HHH(theta_hat) |> MASS::ginv() / 1000
theta_tilde <- theta_hat + Hinv %*% A
theta_tilde
