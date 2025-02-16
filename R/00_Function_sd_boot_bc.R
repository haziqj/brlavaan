##==============================================================================
## Project:     Resampling based bias correction for small sample SEM
##
## Script:      Define function to compute bootstrap bias corrected estimates
##==============================================================================

sd_boot_bc <- function(fit = NULL,
                       ndraws = 500L,
                       boot.resamples = NULL,
                       jack.resamples = NULL){

  # Obtain boot resamples if not provided
  if(is.null(boot.resamples)){
    boot.resamples <- sd_boot_draws(fit, ndraws = ndraws)
  }

  # Extract point estimates per iteration
  u <- boot.resamples[boot.resamples$tryerror != 1,
                      grepl(pattern = "est_", x = colnames(boot.resamples))]

  # Extract SE per iteration
  w <- boot.resamples[boot.resamples$tryerror != 1,
                      grepl(pattern = "se_", x = colnames(boot.resamples))]

  # Compute estimates of standard error, bias, and bias-corrected estimates
  boot.se <- apply(u, MARGIN = 2, stats::sd)
  boot.bias <- colMeans(u) - lavaan::coef(fit)
  theta.bc <- 2*lavaan::coef(fit) - colMeans(u)

  #-----------------------------------------------------------------------------
  # STANDARD NORMAL CONFIDENCE INTERVAL
  # ----------------------------------------------------------------------------

  # Usual formula for calculating CIs, using the normal distribution value 1.96
  # as the multiplier of the standard error. However, there are two differences:
  # --- (a) use our bootstrap estimate of the standard error in the formula
  # --- (b) make an adjustment for the estimated bias

  ci.norm <- c("LB" = lavaan::coef(fit) - boot.bias - 1.96*boot.se,
               "UB" = lavaan::coef(fit) - boot.bias + 1.96*boot.se)

  #-----------------------------------------------------------------------------
  # BASIC BOOTSTRAP CONFIDENCE INTERVAL
  #-----------------------------------------------------------------------------

  # Also known as reverse percentile interval
  # LB = 2*\hat{\theta} - Q_{1-\alpha/2}
  # UB = 2*\hat{\theta} - Q_{\alpha/2}
  # where Q_xxx is the quantile of the bootstrap distribution

  ci.basic <- c("LB" = 2*lavaan::coef(fit) - apply(X = u,
                                           MARGIN = 2,
                                           FUN = stats::quantile,
                                           probs = 0.975,
                                           type = 6),
                "UB" = 2*lavaan::coef(fit) - apply(X = u,
                                           MARGIN = 2,
                                           FUN = stats::quantile,
                                           probs = 0.025,
                                           type = 6))

  #-----------------------------------------------------------------------------
  # Studentized bootstrap confidence interval
  #-----------------------------------------------------------------------------

  # a) compute T_boot = (\hat{\theta}_boot - \hat{\theta})/s_boot
  # b) order all T_boot values from lowest to highest
  # c) choose the 2.5 and 97.5 percentiles as critical values
  # d) LB = \hat{\theta} - T_boot^0.975*SE(\hat{\theta})_boot
  #    UB = \hat{\theta} - T_boot^0.025*SE(\hat{\theta})_boot

  boot.t.table <- sweep(x = u, MARGIN = 2, STATS = lavaan::coef(fit), FUN = "-")/w
  critical.values <- apply(boot.t.table,
                           MARGIN = 2,
                           stats::quantile,
                           probs = c(.025,.975),
                           type = 6,
                           na.rm = TRUE)

  lav.se <-  lavaan::parameterEstimates(fit)[!is.na(lavaan::parameterEstimates(fit)$z), "se"]
  ci.stud <- c("LB" = lavaan::coef(fit) - critical.values[2, ]*lav.se,
               "UB" = lavaan::coef(fit) - critical.values[1, ]*lav.se)

  #-----------------------------------------------------------------------------
  # PERCENTILE CONFIDENCE INTERVAL
  #-----------------------------------------------------------------------------

  # a) order observed bootstrap estimates from lowest to highest
  # b) LB = value at the 2.5th percentile
  #    UB = value at the 97.5th percentile

  LBs <- apply(X = u,
               MARGIN = 2,
               FUN = stats::quantile,
               probs = 0.025,
               type = 6)

  UBs <- apply(X = u,
               MARGIN = 2,
               FUN = stats::quantile,
               probs = 0.975,
               type = 6)

  names(LBs) <- names(UBs) <- names(lavaan::coef(fit))

  ci.perc <- c("LB" = LBs,
               "UB" = UBs)

  #-----------------------------------------------------------------------------
  # BIAS CORRECTED AND ACCELERATED CONFIDENCE INTERVAL
  #-----------------------------------------------------------------------------

  # BCa intervals require estimating a bias term + an acceleration term:
  # b.term = bias estimate (often called z0)
  # a.term = acceleration estimate

  if(is.null(jack.resamples)){
    jack.resamples <- sd_jack_draws(fit)
    jack.est <- jack.resamples[ ,grepl("est_", colnames(jack.resamples))]
  } else {
    jack.est <- jack.resamples[ ,grepl("est_", colnames(jack.resamples))]
  }

  counts <- colMeans(sweep(u, MARGIN = 2, STATS = lavaan::coef(fit), FUN = "<"))
  counts <- ifelse(counts == 0, counts + 1e-09, counts)
  counts <- ifelse(counts == 1, counts - 1e-09, counts)
  b.term <- stats::qnorm(counts)

  diffs <- sweep(-jack.est, MARGIN = 2, STATS = colMeans(jack.est), FUN = "+")

  a.top <- colSums(diffs^3)
  a.bottom <- 6*(colSums(diffs^2)^(3/2))
  a.term <- a.top/a.bottom

  alpha <- 0.025
  alpha1 <- stats::pnorm(b.term + (b.term + stats::qnorm(alpha))/
                    (1-a.term*(b.term+stats::qnorm(alpha))))
  alpha2 <- stats::pnorm(b.term + (b.term + stats::qnorm(1-alpha))/
                    (1-a.term*(b.term + stats::qnorm(1-alpha))))

  LBs <- diag(apply(u,
                    FUN = stats::quantile,
                    MARGIN = 2,
                    probs = alpha1,
                    type = 6))

  UBs <- diag(apply(u,
                    FUN = stats::quantile,
                    MARGIN = 2,
                    probs = alpha2,
                    type = 6))

  names(LBs) <- names(UBs) <- names(lavaan::coef(fit))

  ci.bca <- c("LB" = LBs,
              "UB" = UBs)

  #-----------------------------------------------------------------------------
  # Return output as a list
  #-----------------------------------------------------------------------------

  return(list("boot.resamples" = boot.resamples,
              "boot.se" = boot.se,
              "boot.bias" = boot.bias,
              "theta.bc" = theta.bc,
              "ci.norm" = ci.norm,
              "ci.basic" = ci.basic,
              "ci.stud" = ci.stud,
              "ci.perc" = ci.perc,
              "ci.bca" = ci.bca))

}
