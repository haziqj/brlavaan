library(tidyverse)
library(lavaan)
theme_set(theme_bw())

# Generate two-factor SEM data -------------------------------------------------
n <- 50
dat <- simulateData(
  model = "
    eta1 =~ 1*y1 + 0.8*y2 + 0.6*y3
    eta2 =~ 1*y4 + 0.8*y5 + 0.6*y6
    eta2 ~ 0.3*eta1
  ",
  sample.nobs = n
)

theta_true <- c(
  0.8, 0.6, 0.8, 0.6,  # lambda
  0.3,  # beta
  rep(1, 6),  # theta_eps
  rep(1, 2)  # psi
)
m <- length(theta_true)

# {lavaan} fit -----------------------------------------------------------------
mod <- "
    eta1 =~ y1 + y2 + y3
    eta2 =~ y4 + y5 + y6
    eta2 ~ eta1
  "
fit_lav <- sem(mod, dat)
theta_lav <- coef(fit_lav)

# ggplot(data.frame(x = theta_lav, y = theta_true), aes(x, y)) +
#   geom_point() +
#   geom_abline()

# Maximum likelihood (manual ish) ----------------------------------------------
loglik <- function(theta, whichlog = 6:m) {

  theta[whichlog] <- exp(theta[whichlog])

  suppressWarnings({
    fit <- sem(
      model = mod,
      data = dat,
      start = theta,
      optim.method = "none",
      baseline = FALSE
    )
  as.numeric(logLik(fit))
  })
}

# theta_ml <- optim(
#   par = rnorm(m),
#   fn = loglik,
#   method = "BFGS",
#   control = list(fnscale = -1)
# )$par
# theta_ml[6:m] <- exp(theta_ml[6:m])
#
# sum((theta_ml - theta_lav) ^ 2)
# [1] 7.99439e-09

# Implicit bias reduction ------------------------------------------------------
JJJ <- function(theta, object = fit_lav, whichlog = 6:m) {
  # Note: This is the e(theta) matrix

  theta[whichlog] <- exp(theta[whichlog])

  # extract slots
  lavmodel       <- object@Model
  lavsamplestats <- object@SampleStats
  lavdata        <- object@Data
  lavoptions     <- object@Options

  # update lavmodel with 'new' set of values
  lavmodel <- lav_model_set_parameters(lavmodel, x = theta)

  lavaan:::lav_model_information_firstorder(
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
}

HHH <- function(theta, object = fit_lav, whichlog = 6:m,
                unit_information = FALSE) {
  # Note: This is the j(theta) matrix

  theta[whichlog] <- exp(theta[whichlog])

  # extract slots
  lavmodel       <- object@Model
  lavsamplestats <- object@SampleStats
  lavdata        <- object@Data
  lavoptions     <- object@Options
  lavcache       <- object@Cache
  lavimplied     <- object@implied
  lavh1          <- object@h1

  # update lavmodel with 'new' set of values
  lavmodel <- lav_model_set_parameters(lavmodel, x = theta)

  out <- lavaan:::lav_model_hessian(
    lavmodel       = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata        = lavdata,
    lavoptions     = lavoptions,
    lavcache       = lavcache,
    ceq.simple     = FALSE
  )

  if (isTRUE(unit_information)) out <- out / lavsamplestats@ntotal
  out
}

obj_fun <- function(theta) {
  e <- JJJ(theta)
  j <- HHH(theta)
  # out <- loglik(theta) - 0.5 * sum(diag(MASS::ginv(j) %*% e))
  out <- loglik(theta) - 0.5 * sum(diag(solve(j, e))) - sum(theta ^ 2) / n
  out
}

theta_start <- coef(fit_lav)
theta_start[6:m] <- log(theta_start[6:m])

theta_iRBMp <- optim(
  par = theta_start,
  fn = obj_fun,
  method = "BFGS",
  # lower = c(rep(-Inf, 5), rep(-2, m - 5)),
  # upper = c(rep(Inf, 5), rep(2, m - 5)),
  control = list(fnscale = -1, trace = 1)
)$par
theta_iRBMp[6:m] <- exp(theta_iRBMp[6:m])

ggplot(data.frame(x = theta_lav, y = theta_iRBMp), aes(x, y))+
  geom_point() +
  geom_abline()

