##==============================================================================
## Project:     Resampling based bias correction for small sample SEM
##
## Script:      Define function to compute jackknife bias corrected estimates
##==============================================================================

sd_jack_bc <- function(fit = NULL, jack.resamples = NULL) {
  
  # Obtain boot resamples if not provided
  if(is.null(jack.resamples)){
    jack.resamples <- sd_jack_draws(fit)
  }
  
  # Determine sample size
  n <- fit@Data@nobs[[1]]
  
  # Compute jackknife estimates
  if(sum(jack.resamples[ , "convergence"], na.rm = TRUE) == n) {
    
    u <- jack.resamples[ , paste0("est_", names(coef(fit)))]
    jack.se <- sqrt((n-1)/n*colSums(sweep(u, 2, colMeans(u), FUN = "-")^2))
    jack.bias <- (n-1)*(colMeans(u)-coef(fit))
    theta.bc <- n*coef(fit) - (n-1)*colMeans(u)
    
    ci.norm <- c("LB" = coef(fit) - jack.bias - 1.96*jack.se,
                 "UB" = coef(fit) - jack.bias + 1.96*jack.se)

  } else {
    
    jack.se <- rep(NA, times = length(coef(fit)))
    jack.bias <- rep(NA, times = length(coef(fit)))
    theta.bc <- rep(NA, times = length(coef(fit)))
    ci.norm <- rep(NA, times = length(coef(fit))*2)
    
  }
  
  # Return values as list
  return(list("jack.resamples" = jack.resamples,
              "jack.se" = jack.se,
              "jack.bias" = jack.bias,
              "theta.bc" = theta.bc, 
              "ci.norm" = ci.norm))
  
}

