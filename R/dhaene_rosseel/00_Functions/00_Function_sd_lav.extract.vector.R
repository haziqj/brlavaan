##==============================================================================
## Project:     Resampling based bias correction for small sample SEM
##
## Script:      Define function to extract all necessary info from a lavaan fit
##==============================================================================

sd_lav.extract.vector <- function(fit = NULL, extras = FALSE){
  
  # extract sample size
  nobs <- fit@Data@nobs[[1]]
  
  # count how many unique cases
  unique.obs <- nrow(unique(fit@Data@X[[1]]))
  
  # check if cov matrix is positive definite
  cov.posdef <- matrixcalc::is.positive.definite(fit@SampleStats@cov[[1]])
  chol2inv.fail <- inherits(x = try(chol2inv(chol(fit@SampleStats@cov[[1]])), 
                                    silent = TRUE), what = "try-error")
  
  # determine for each parameter which bound was hit (if any)
  # 0 = no hits, 1 = lower bound, 2 = upper bound
  PT <- parTable(fit)[parTable(fit)$free != 0, ]
  lb.hits <- ifelse(PT$est == PT$lower, yes = 1, no = 0)
  ub.hits <- ifelse(PT$est == PT$upper, yes = 2, no = 0)
  hit.bound <- replace(x = lb.hits, 
                       list = which(lb.hits == 0), 
                       values = ub.hits[which(lb.hits == 0)])
  
  # extract convergence, coefficients, SE & confidence limits
  PE <- parameterEstimates(fit)[as.numeric(rownames(PT)), ] %>%
    tidyr::unite("param", c("lhs", "op", "rhs"), remove = TRUE, sep = "")
  
  if(extras == TRUE) {
    
    info <- c(lavInspect(fit, what = "converged"), 
              nobs,
              unique.obs,
              cov.posdef, 
              chol2inv.fail,
              PE$est, 
              PE$se,
              PE$ci.lower,
              PE$ci.upper,
              hit.bound)
    
    names(info) <- c("convergence",
                     "nobs", 
                     "unique.obs", 
                     "cov.posdef",
                     "chol2inv.fail",
                     paste0("est_", PE$param), 
                     paste0("se_", PE$param),
                     paste0("ci.lower_", PE$param), 
                     paste0("ci.upper_", PE$param), 
                     paste0("hit.bound_", PE$param))
    
  } else {
    
    info <- c(lavInspect(fit, what = "converged"), 
              nobs,
              PE$est, 
              PE$se,
              PE$ci.lower,
              PE$ci.upper,
              hit.bound)
    
    names(info) <- c("convergence",
                     "nobs",
                     paste0("est_", PE$param), 
                     paste0("se_", PE$param),
                     paste0("ci.lower_", PE$param), 
                     paste0("ci.upper_", PE$param), 
                     paste0("hit.bound_", PE$param))
    
  }
 
  return (info)
  
}
