imp_twofac <- function(dat, trace = FALSE) {
  # {lavaan} fit ---------------------------------------------------------------
  mod <- "
    eta1 =~ y1 + y2 + y3
    eta2 =~ y4 + y5 + y6
  "
  fit_lav <- sem(mod, dat)
  theta_lav <- coef(fit_lav)

  # Functions ------------------------------------------------------------------
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
    # print(theta)
    e <- JJJ(theta)
    j <- HHH(theta) + diag(1e-7, nrow(e))
    # out <- loglik(theta) - 0.5 * sum(diag(MASS::ginv(j) %*% e))
    out <- loglik(theta) - 0.5 * sum(diag(solve(j, e))) - sum(theta ^ 2) / n
    out
  }

  # Start iRBMp fit ------------------------------------------------------------
  theta_start <- theta_lav
  theta_start[6:m] <- log(abs(theta_start[6:m]))

  res <- try(optim(
    par = theta_start,
    fn = obj_fun,
    method = "L-BFGS-B",
    lower = rep(-5, m),
    upper = rep(5, m),
    control = list(fnscale = -1, trace = trace)
  ))
  if (inherits(res, "try-error")) {
    theta_iRBMp <- res
  } else {
    theta_iRBMp <- res$par
    theta_iRBMp[6:m] <- exp(theta_iRBMp[6:m])
  }

  list(
    ml = theta_lav,
    iRBMp = theta_iRBMp
  )
}
