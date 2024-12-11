##==============================================================================
## Project:     Resampling based bias correction for small sample SEM
##
## Script:      Define function to extract REML estimates
##==============================================================================

#-------------------------------------------------------------------------------
# Version 1: using merDeriv package
#-------------------------------------------------------------------------------

est.REML <- function(lme4model) {
  
  vc <- as.data.frame(VarCorr(lme4model))$vcov
  fe <- fixef(lme4model)
  
  est <- c("i~~i" = vc[1],
           "i~1" = fe[[1]],
           "s~~s" = vc[2],
           "s~1" = fe[[2]],
           "i~~s" = vc[3],
           "Day1~~Day1" = vc[4],
           "Day2~~Day2" = vc[4],
           "Day3~~Day3" = vc[4],
           "Day4~~Day4" = vc[4],
           "Day5~~Day5" = vc[4],
           "Day6~~Day6" = vc[4],
           "Day7~~Day7" = vc[4],
           "Day8~~Day8" = vc[4],
           "Day9~~Day9" = vc[4],
           "Day10~~Day10" = vc[4])
  
  ses <- sqrt(diag(vcov(lme4model, full = TRUE)))
  
  se <- c("i~~i" = ses[3],
          "i~1" = ses[1],
          "s~~s" = ses[5],
          "s~1" = ses[2],
          "i~~s" = ses[4],
          "Day1~~Day1" = ses[6],
          "Day2~~Day2" = ses[6],
          "Day3~~Day3" = ses[6],
          "Day4~~Day4" = ses[6],
          "Day5~~Day5" = ses[6],
          "Day6~~Day6" = ses[6],
          "Day7~~Day7" = ses[6],
          "Day8~~Day8" = ses[6],
          "Day9~~Day9" = ses[6],
          "Day10~~Day1°" = ses[6])
  
  return(list("est" = est,
              "se" = se))
  
}

#-------------------------------------------------------------------------------
# Version 2: output lme4: only SEs for fixed effects
#-------------------------------------------------------------------------------

est.REML.v2 <- function(lme4model) {
  
  vc <- as.data.frame(VarCorr(lme4model))$vcov
  fe <- fixef(lme4model)
  
  est <- c("i~~i" = vc[1],
           "i~1" = fe[[1]],
           "s~~s" = vc[2],
           "s~1" = fe[[2]],
           "i~~s" = vc[3],
           "Day1~~Day1" = vc[4],
           "Day2~~Day2" = vc[4],
           "Day3~~Day3" = vc[4],
           "Day4~~Day4" = vc[4],
           "Day5~~Day5" = vc[4],
           "Day6~~Day6" = vc[4],
           "Day7~~Day7" = vc[4],
           "Day8~~Day8" = vc[4],
           "Day9~~Day9" = vc[4],
           "Day10~~Day10" = vc[4])
  
  out <- summary(fit.REML)
  ses <- out$coefficients[ , 2]

  se <- c("i~~i" = NA,
          "i~1" = ses[1],
          "s~~s" = NA,
          "s~1" = ses[2],
          "i~~s" = NA,
          "Day1~~Day1" = NA,
          "Day2~~Day2" = NA,
          "Day3~~Day3" = NA,
          "Day4~~Day4" = NA,
          "Day5~~Day5" = NA,
          "Day6~~Day6" = NA,
          "Day7~~Day7" = NA,
          "Day8~~Day8" = NA,
          "Day9~~Day9" = NA,
          "Day10~~Day1°" = NA)
  
  return(list("est" = est,
              "se" = se))
  
}