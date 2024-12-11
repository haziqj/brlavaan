##==============================================================================
## Project:     Resampling based bias correction for small sample SEM
##
## Script:      Define function to perform a simple rescaling on the variance 
##              components (e.g., N/(N-1))
##==============================================================================

sd_varcomp_rescaling <- function(fit = NULL, bc.factor = 1) {
  
  n <- fit@Data@nobs[[1]]
  bc <- n/(n-bc.factor)
  
  PT.bc <- parameterTable(fit) %>%
    filter(free != 0) %>%
    mutate(est.bc = ifelse(op == "~~", est * bc, est))
  
  return(PT.bc$est.bc)
  
}
