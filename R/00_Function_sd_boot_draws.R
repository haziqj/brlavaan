##==============================================================================
## Project:     Resampling based bias correction for small sample SEM
##
## Script:      Define function to resample from original sample and fit the
##              model to the resampled data
##==============================================================================

sd_boot_draws <- function(fit, ndraws = 500L) {

  # extract model options
  opt.means <- fit@Options$meanstructure
  opt.model <- fit@Options$model.type
  opt.bounds <- fit@Options$bounds

  # original sample
  data <- fit@Data@X[[1]]
  colnames(data) <- fit@Data@ov.names[[1]]

  # number of observations
  nobs <- fit@Data@nobs[[1]]

  # save all info resamples
  param.names <- names(lavaan::coef(fit))
  info <- c("resample",
            "tryerror",
            "convergence",
            "nobs",
            "unique.obs",
            "cov.posdef",
            "chol2inv.fail",
            paste0("est_", param.names),
            paste0("se_", param.names),
            paste0("ci.lower_", param.names),
            paste0("ci.upper_", param.names),
            # paste0("z_", param.names),
            # paste0("pvalue_", param.names),
            paste0("hit.bound_", param.names))

  est <- as.data.frame(matrix(NA, nrow = ndraws, ncol = length(info)))
  colnames(est) <- info

  # compute bootstrap estimates
  for(i in seq_len(ndraws)) {

    # draw bootstrap sample
    seed.boot <- 123 + i
    set.seed(seed.boot)
    Y <- data[sample.int(nrow(data), size = nrow(data), replace = TRUE), ]

    # fit model
    if(opt.model == "growth"){
      fit.i <- try(lavaan::growth(model = lavaan::parameterTable(fit),
                          data = Y,
                          estimator = "ML",
                          optim.attempts = 1,
                          bounds = opt.bounds,
                          meanstructure = opt.means), silent = TRUE)

    } else {
      fit.i <- try(lavaan::sem(model = lavaan::parameterTable(fit),
                       data = Y,
                       estimator = "ML",
                       optim.attempts = 1,
                       bounds = opt.bounds,
                       meanstructure = opt.means), silent = TRUE)

    }

    if(!inherits(fit.i, "try-error")){

      # register try-error
      tryerror <- 0

      # extract all necessary info
      all.info <- sd_lav.extract.vector(fit.i, extras = TRUE)

      # fill est matrix
      est[i, ] <- c("resample" = i,
                    "tryerror" = tryerror,
                    all.info)

    } else {

      # register try-error
      tryerror <- 1
      cov.matrix <- stats::cov(Y)*(nobs-1)/nobs
      cov.posdef <- matrixcalc::is.positive.definite(cov.matrix)
      chol2inv.fail <- inherits(x = try(chol2inv(chol(cov.matrix)),
                                        silent = TRUE), what = "try-error")

      est[i, ] <- c("resample" = i,
                    "tryerror" = tryerror,
                    "convergence" = 0,
                    "nobs" = nrow(Y),
                    "unique.nobs" = nrow(unique(Y)),
                    "cov.posdef" = cov.posdef,
                    "chol2inv.fail" = chol2inv.fail,
                    rep(NA, times = ncol(est) - 7))
    }

  }

  return("boot.resamples" = est)

}
