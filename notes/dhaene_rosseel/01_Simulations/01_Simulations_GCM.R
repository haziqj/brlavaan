##==============================================================================
## Project:     Resampling based bias correction for small sample SEM
##
## Script:      Simulations growth curve model
##==============================================================================

#-------------------------------------------------------------------------------
# Set parameters for this script
#-------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
(nobs <- as.numeric(args[1]))
(rel <- args[2])
(dist <- args[3])
(nsim <- as.numeric(args[4]))
(jobid <- as.numeric(args[5]))
ndraws <- 500

# Specify estimation methods
est.methods <- c("MLB" = 1, 
                 "JB" = 2, 
                 "BB" = 3,
                 "Ozenne" = 4,
                 "simple.min1" = 5,
                 "simple.min2" = 6,
                 "simple.min3" = 7, 
                 "REML" = 8)

conf.methods <- c("MLB" = 1, 
                  "JB" = 2, 
                  "BB-ci.norm" = 3, 
                  "BB-ci.basic" = 4,
                  "BB-ci.stud" = 5,
                  "BB-ci.perc" = 6,
                  "BB-ci.bca" = 7,
                  "Ozenne" = 8,
                  "simple.min1" = 9,
                  "simple.min2" = 10,
                  "simple.min3" = 11,
                  "REML" = 12)

#-------------------------------------------------------------------------------
# Load libraries
#-------------------------------------------------------------------------------

library(lavaan)
library(dplyr)
library(lme4)
library(merDeriv)

#-------------------------------------------------------------------------------
# Source files
#-------------------------------------------------------------------------------

source("00_Functions/00_Function_sd_lav.extract.vector.R")
source("00_Functions/00_Function_sd_jack_draws.R")
source("00_Functions/00_Function_sd_boot_draws.R")
source("00_Functions/00_Function_sd_jack_bc.R")
source("00_Functions/00_Function_sd_boot_bc.R")
source("00_Functions/00_Function_oz_bc_growth.R")
source("00_Functions/00_Function_sd_varcomp_rescaling.R")
source("00_Functions/00_Function_REML.R")

#-------------------------------------------------------------------------------
# Set up data generating model
#-------------------------------------------------------------------------------

if(readr::parse_number(rel)/100 == 0.8) {
  
  # Population model in matrices
  n.factors <- 2
  n.indicators <- 10
  lambda <- matrix(data = c(rep(1, times = n.indicators), 0:(n.indicators-1)),
                   nrow = n.indicators, 
                   ncol = n.factors)
  psi <- matrix(data = c(550, 40, 40, 100), 
                nrow = n.factors, 
                ncol = n.factors)
  theta <- matrix(data = 0, 
                  nrow = n.indicators,
                  ncol = n.indicators)
  diag(theta) <- 500
  
  # Population model in lavaan syntax
  pop.model <-    '   # intercept with coefficients fixed to 1
                      i =~  1*Day0 + 1*Day1 + 1*Day2 + 1*Day3 + 1*Day4 +
                            1*Day5 + 1*Day6 + 1*Day7 + 1*Day8 + 1*Day9
                                
                      # slope with coefficients fixed to 0:9 (number of days)
                      s =~  0*Day0 + 1*Day1 + 2*Day2 + 3*Day3 + 4*Day4 +
                            5*Day5 + 6*Day6 + 7*Day7 + 8*Day8 + 9*Day9
  
                      i ~~ 550*i
                      i ~ 0*1
                      
                      s ~~ 100*s
                      s ~ 0*1
                      
                      i ~~ 40*s
                      
                      Day0 ~~ 500*Day0
                      Day1 ~~ 500*Day1
                      Day2 ~~ 500*Day2
                      Day3 ~~ 500*Day3
                      Day4 ~~ 500*Day4
                      Day5 ~~ 500*Day5
                      Day6 ~~ 500*Day6
                      Day7 ~~ 500*Day7
                      Day8 ~~ 500*Day8
                      Day9 ~~ 500*Day9 '
  
  # Analysis model in lavaan syntax
  model <-        '   # intercept with coefficients fixed to 1
                      i =~  1*Day0 + 1*Day1 + 1*Day2 + 1*Day3 + 1*Day4 +
                            1*Day5 + 1*Day6 + 1*Day7 + 1*Day8 + 1*Day9
                                
                      # slope with coefficients fixed to 0:9 (number of days)
                      s =~  0*Day0 + 1*Day1 + 2*Day2 + 3*Day3 + 4*Day4 +
                            5*Day5 + 6*Day6 + 7*Day7 + 8*Day8 + 9*Day9
  
                      i ~~ start(550)*i
                      i ~ start(0)*1
                      
                      s ~~ start(100)*s
                      s ~ start(0)*1
                      
                      i ~~ start(40)*s
                      
                      # apply equality constraints
                      Day0 ~~ v*Day0 + start(500)*Day0
                      Day1 ~~ v*Day1 + start(500)*Day1
                      Day2 ~~ v*Day2 + start(500)*Day2
                      Day3 ~~ v*Day3 + start(500)*Day3
                      Day4 ~~ v*Day4 + start(500)*Day4
                      Day5 ~~ v*Day5 + start(500)*Day5
                      Day6 ~~ v*Day6 + start(500)*Day6
                      Day7 ~~ v*Day7 + start(500)*Day7
                      Day8 ~~ v*Day8 + start(500)*Day8
                      Day9 ~~ v*Day9 + start(500)*Day9 '

} else if(readr::parse_number(rel)/100 == 0.5) {
  
  # Population model in matrices
  n.factors <- 2
  n.indicators <- 10
  lambda <- matrix(data = c(rep(1, times = n.indicators), 0:(n.indicators-1)),
                   nrow = n.indicators, 
                   ncol = n.factors)
  psi <- matrix(data = c(275, 20, 20, 50), 
                nrow = n.factors, 
                ncol = n.factors)
  theta <- matrix(data = 0, 
                  nrow = n.indicators,
                  ncol = n.indicators)
  diag(theta) <- 1300
  
  # Population model in lavaan syntax
  pop.model <-    '   # intercept with coefficients fixed to 1
                    i =~  1*Day0 + 1*Day1 + 1*Day2 + 1*Day3 + 1*Day4 +
                          1*Day5 + 1*Day6 + 1*Day7 + 1*Day8 + 1*Day9
                              
                    # slope with coefficients fixed to 0:9 (number of days)
                    s =~  0*Day0 + 1*Day1 + 2*Day2 + 3*Day3 + 4*Day4 +
                          5*Day5 + 6*Day6 + 7*Day7 + 8*Day8 + 9*Day9

                    i ~~ 275*i
                    i ~ 0*1
                    
                    s ~~ 50*s
                    s ~ 0*1
                    
                    i ~~ 20*s
                    
                    Day0 ~~ 1300*Day0
                    Day1 ~~ 1300*Day1
                    Day2 ~~ 1300*Day2
                    Day3 ~~ 1300*Day3
                    Day4 ~~ 1300*Day4
                    Day5 ~~ 1300*Day5
                    Day6 ~~ 1300*Day6
                    Day7 ~~ 1300*Day7
                    Day8 ~~ 1300*Day8
                    Day9 ~~ 1300*Day9 '
  
  # Analysis model in lavaan syntax
  model <-        '   # intercept with coefficients fixed to 1
                      i =~  1*Day0 + 1*Day1 + 1*Day2 + 1*Day3 + 1*Day4 +
                            1*Day5 + 1*Day6 + 1*Day7 + 1*Day8 + 1*Day9
                                
                      # slope with coefficients fixed to 0:9 (number of days)
                      s =~  0*Day0 + 1*Day1 + 2*Day2 + 3*Day3 + 4*Day4 +
                            5*Day5 + 6*Day6 + 7*Day7 + 8*Day8 + 9*Day9
  
                      i ~~ start(275)*i
                      i ~ start(0)*1
                      
                      s ~~ start(50)*s
                      s ~ start(0)*1
                      
                      i ~~ start(20)*s
                      
                      # apply equality constraints
                      Day0 ~~ v*Day0 + start(1300)*Day0
                      Day1 ~~ v*Day1 + start(1300)*Day1
                      Day2 ~~ v*Day2 + start(1300)*Day2
                      Day3 ~~ v*Day3 + start(1300)*Day3
                      Day4 ~~ v*Day4 + start(1300)*Day4
                      Day5 ~~ v*Day5 + start(1300)*Day5
                      Day6 ~~ v*Day6 + start(1300)*Day6
                      Day7 ~~ v*Day7 + start(1300)*Day7
                      Day8 ~~ v*Day8 + start(1300)*Day8
                      Day9 ~~ v*Day9 + start(1300)*Day9 '
  
} else {cat("'rel' needs to be 0.5 or 0.8")}

#-------------------------------------------------------------------------------
# Prepare containers to store estimates
#-------------------------------------------------------------------------------

# Test model fit (shortcut to extract all info)
fit <- growth(data = simulateData(pop.model, empirical = TRUE, sample.nobs = 15), 
              model = model, 
              meanstructure = TRUE, 
              optim.attempts = 1,
              bounds = "standard")

# Define what needs to be stored per container
all.info <- list("estimates" = c("jobid", 
                                 "iteration", 
                                 "seed", 
                                 "method",
                                 "ci.type",
                                 "tryerror",
                                 names(sd_lav.extract.vector(fit))),
                 "resamples" = c("jobid", 
                                 "iteration", 
                                 "seed", 
                                 "resample",
                                 "tryerror",
                                 names(sd_lav.extract.vector(fit, TRUE))))

# Estimates
est <- as.data.frame(matrix(data = NA, 
                            nrow = nsim * length(conf.methods), 
                            ncol = length(all.info$estimates), 
                            dimnames = list(NULL, all.info$estimates)))

# Jackknife resamples 
jack.resamples <- as.data.frame(matrix(data = NA, 
                                       nrow = nsim*nobs, 
                                       ncol = length(all.info$resamples), 
                                       dimnames = list(NULL, all.info$resamples)))

# Bootstrap resamples 
boot.resamples <- as.data.frame(matrix(data = NA, 
                                       nrow = nsim*ndraws, 
                                       ncol = length(all.info$resamples), 
                                       dimnames = list(NULL, all.info$resamples)))

#-------------------------------------------------------------------------------
# Actual simulation loop
#-------------------------------------------------------------------------------

(starting <- Sys.time())
cat("Ready to loop over ", nsim, " simulations! \n")

for(i in 1:nsim) {
  
  #-----------------------------------------------------------------------------
  # Data generation
  #-----------------------------------------------------------------------------
  
  seed <- 1234 + nsim*(jobid-1) + i
  set.seed(seed)
  
  if(dist == "Normal") {
    
    factors <- t(covsim::rIG(N = nobs, 
                           sigma.target = psi, 
                           skewness = rep(0, times = ncol(lambda)),
                           excesskurtosis = rep(0, times = ncol(lambda)),
                           reps = 1)[[1]])
    
    e <- t(covsim::rIG(N = nobs, 
                       sigma.target = theta, 
                       skewness = rep(0, times = nrow(theta)), 
                       excesskurtosis = rep(0, times = nrow(theta)), 
                       reps = 1)[[1]])
    
  } else if(dist == "Kurtosis") {
    
    factors <- t(covsim::rIG(N = nobs, 
                             sigma.target = psi, 
                             skewness = rep(0, times = ncol(lambda)),
                             excesskurtosis = rep(6, times = ncol(lambda)),
                             reps = 1)[[1]])
    
    e <- t(covsim::rIG(N = nobs, 
                       sigma.target = theta, 
                       skewness = rep(0, times = nrow(theta)), 
                       excesskurtosis = rep(6, times = nrow(theta)), 
                       reps = 1)[[1]])
    
  } else if(dist == "NonNormal") {
    
    factors <- t(covsim::rIG(N = nobs, 
                             sigma.target = psi, 
                             skewness = rep(-2, times = ncol(lambda)),
                             excesskurtosis = rep(6, times = ncol(lambda)),
                             reps = 1)[[1]])
    
    e <- t(covsim::rIG(N = nobs, 
                       sigma.target = theta, 
                       skewness = rep(-2, times = nrow(theta)), 
                       excesskurtosis = rep(6, times = nrow(theta)), 
                       reps = 1)[[1]])
    
    
  } else {cat("'dist' needs to be 'Normal', 'Kurtosis' or 'NonNormal'")}
  
  data <- as.data.frame(t(lambda %*% factors + e))
  colnames(data) <- paste0("Day", 0:9)
  
  #-----------------------------------------------------------------------------
  # Fit model
  #-----------------------------------------------------------------------------
  
  fit <- try(growth(data = data, 
                    model = model, 
                    meanstructure = TRUE, 
                    optim.attempts = 1,
                    bounds = "standard"), silent = TRUE)
  
  if(!inherits(fit, "try-error")){
    
    # Maximum Likelihood Bounds (MLB)
    est[nsim*0 + i, ] <- c("jobid" = jobid, 
                           "iteration" = i,
                           "seed" = seed,
                           "method" = est.methods["MLB"],
                           "ci.type" = conf.methods["MLB"],
                           "tryerror" = 0,
                           sd_lav.extract.vector(fit))
    
    cat("Finished: i = ", i, "\t --- MLB", "\n", sep = "")
    
    # Jackknife Bounds: bias-corrected estimates (JB)
    jack.out <- sd_jack_bc(fit = fit)
    
    a <- (i-1)*nobs + 1
    b <- i*nobs
    jack.resamples[a:b, ] <- cbind("jobid" = jobid, 
                                   "iteration" = i, 
                                   "seed" = seed,
                                   jack.out$jack.resamples)
    
    est[nsim*1 + i, ] <- c("jobid" = jobid, 
                           "iteration" = i,
                           "seed" = seed,
                           "method" = est.methods["JB"],
                           "ci.type" = conf.methods["JB"],
                           "tryerror" = 0,
                           "convergence" = lavInspect(fit, "converged"),
                           "nobs" = nobs,
                           jack.out$theta.bc,
                           jack.out$jack.se,
                           jack.out$ci.norm,
                           rep(NA, times = length(coef(fit))))  
    
    cat("Finished: i = ", i, "\t --- JB", "\n", sep = "")
    
    # Bootstrap Bounds: bias-corrected estimates (BB) - ci.norm
    boot.out <- sd_boot_bc(fit = fit, 
                           ndraws = ndraws,
                           jack.resamples = jack.out$jack.resamples)
    
    a <- (i-1)*ndraws + 1
    b <- i*ndraws
    boot.resamples[a:b, ] <- cbind("jobid" = jobid, 
                                   "iteration" = i, 
                                   "seed" = seed,
                                   boot.out$boot.resamples)
    
    est[nsim*2 + i, ] <- c("jobid" = jobid, 
                           "iteration" = i,
                           "seed" = seed,
                           "method" = est.methods["BB"],
                           "ci.type" = conf.methods["BB-ci.norm"],
                           "tryerror" = 0,
                           "convergence" = lavInspect(fit, "converged"),
                           "nobs" = nobs,
                           boot.out$theta.bc,
                           boot.out$boot.se,
                           boot.out$ci.norm,
                           rep(NA, times = length(coef(fit))))  
    
    # Bootstrap Bounds: bias-corrected estimates (BB) - ci.basic
    est[nsim*3 + i, ] <- c("jobid" = jobid, 
                           "iteration" = i,
                           "seed" = seed,
                           "method" = est.methods["BB"],
                           "ci.type" = conf.methods["BB-ci.basic"],
                           "tryerror" = 0,
                           "convergence" = lavInspect(fit, "converged"),
                           "nobs" = nobs,
                           boot.out$theta.bc,
                           boot.out$boot.se,
                           boot.out$ci.basic,
                           rep(NA, times = length(coef(fit))))  
    
    # Bootstrap Bounds: bias-corrected estimates (BB) - ci.stud
    est[nsim*4 + i, ] <- c("jobid" = jobid, 
                           "iteration" = i,
                           "seed" = seed,
                           "method" = est.methods["BB"],
                           "ci.type" = conf.methods["BB-ci.stud"],
                           "tryerror" = 0,
                           "convergence" = lavInspect(fit, "converged"),
                           "nobs" = nobs,
                           boot.out$theta.bc,
                           boot.out$boot.se,
                           boot.out$ci.stud,
                           rep(NA, times = length(coef(fit))))  
    
    # Bootstrap Bounds: bias-corrected estimates (BB) - ci.perc
    est[nsim*5 + i, ] <- c("jobid" = jobid, 
                           "iteration" = i,
                           "seed" = seed,
                           "method" = est.methods["BB"],
                           "ci.type" = conf.methods["BB-ci.perc"],
                           "tryerror" = 0,
                           "convergence" = lavInspect(fit, "converged"),
                           "nobs" = nobs,
                           boot.out$theta.bc,
                           boot.out$boot.se,
                           boot.out$ci.perc,
                           rep(NA, times = length(coef(fit))))  
    
    # Bootstrap Bounds: bias-corrected estimates (BB) - ci.bca
    est[nsim*6 + i, ] <- c("jobid" = jobid, 
                           "iteration" = i,
                           "seed" = seed,
                           "method" = est.methods["BB"],
                           "ci.type" = conf.methods["BB-ci.bca"],
                           "tryerror" = 0,
                           "convergence" = lavInspect(fit, "converged"),
                           "nobs" = nobs,
                           boot.out$theta.bc,
                           boot.out$boot.se,
                           boot.out$ci.bca,
                           rep(NA, times = length(coef(fit))))  
    
    cat("Finished: i = ", i, "\t --- BB", "\n", sep = "")
    
    # Ozenne bias-corrected estimates (Ozenne)
    ozenne.out <- try(yr_ozenne_bcs(fit, return.se = TRUE), silent = FALSE)
    
    if(!inherits(ozenne.out, "try-error")){
      
      ci.norm.oz <- c("LB" = ozenne.out$est.corrected - 1.96*ozenne.out$se,
                      "UB" = ozenne.out$est.corrected + 1.96*ozenne.out$se)
      
      est[nsim*7 + i, ] <- c("jobid" = jobid, 
                             "iteration" = i,
                             "seed" = seed,
                             "method" = est.methods["Ozenne"],
                             "ci.type" = conf.methods["Ozenne"],
                             "tryerror" = 0,
                             "convergence" = lavInspect(fit, "converged"),
                             "nobs" = nobs,
                             ozenne.out$est.corrected,
                             ozenne.out$se,
                             ci.norm.oz,
                             rep(NA, times = length(coef(fit))))  
      
    } else {
      
      est[nsim*7 + i, ] <- c("jobid" = jobid, 
                             "iteration" = i,
                             "seed" = seed,
                             "method" = est.methods["Ozenne"],
                             "ci.type" = conf.methods["Ozenne"],
                             "tryerror" = 1,
                             "convergence" = NA,
                             "nobs" = nobs,
                             rep(NA, times = length(all.info$estimates)-8))    
    }
    
    cat("Finished: i = ", i, "\t --- Ozenne", "\n", sep = "")
    
    # Simple bias-corrected estimates (simple): n/n-1
    est_simple.bc <- sd_varcomp_rescaling(fit, bc.factor = 1)
    lav.se <- sd_lav.extract.vector(fit)[grep("se_", 
                                              names(sd_lav.extract.vector(fit)))]
    
    ci.norm.simple <- c("LB" = est_simple.bc - 1.96*lav.se,
                        "UB" = est_simple.bc + 1.96*lav.se)
    
    est[nsim*8 + i, ] <- c("jobid" = jobid, 
                           "iteration" = i,
                           "seed" = seed,
                           "method" = est.methods["simple.min1"],
                           "ci.type" = conf.methods["simple.min1"],
                           "tryerror" = 0,
                           "convergence" = lavInspect(fit, "converged"),
                           "nobs" = nobs,
                           est_simple.bc,
                           lav.se,
                           ci.norm.simple,
                           rep(NA, times = length(coef(fit))))
    
    cat("Finished: i = ", i, "\t --- simple n/n-1", "\n", sep = "")
    
    # Simple bias-corrected estimates (simple): n/n-2
    est_simple.bc.2 <- sd_varcomp_rescaling(fit, bc.factor = 2)
    ci.norm.simple.2 <- c("LB" = est_simple.bc.2 - 1.96*lav.se,
                          "UB" = est_simple.bc.2 + 1.96*lav.se)
    
    est[nsim*9 + i, ] <- c("jobid" = jobid, 
                           "iteration" = i,
                           "seed" = seed,
                           "method" = est.methods["simple.min2"],
                           "ci.type" = conf.methods["simple.min2"],
                           "tryerror" = 0,
                           "convergence" = lavInspect(fit, "converged"),
                           "nobs" = nobs,
                           est_simple.bc.2,
                           lav.se,
                           ci.norm.simple.2,
                           rep(NA, times = length(coef(fit))))  
    
    cat("Finished: i = ", i, "\t --- simple n/n-2", "\n", sep = "")
    
    # Simple bias-corrected estimates (simple): n/n-3
    est_simple.bc.3 <- sd_varcomp_rescaling(fit, bc.factor = 3)
    ci.norm.simple.3 <- c("LB" = est_simple.bc.3 - 1.96*lav.se,
                          "UB" = est_simple.bc.3 + 1.96*lav.se)
    
    est[nsim*10 + i, ] <- c("jobid" = jobid, 
                            "iteration" = i,
                            "seed" = seed,
                            "method" = est.methods["simple.min3"],
                            "ci.type" = conf.methods["simple.min3"],
                            "tryerror" = 0,
                            "convergence" = lavInspect(fit, "converged"),
                            "nobs" = nobs,
                            est_simple.bc.3,
                            lav.se,
                            ci.norm.simple.3,
                            rep(NA, times = length(coef(fit))))  
    
    cat("Finished: i = ", i, "\t --- simple n/n-3", "\n", sep = "")
    
    # REML
    data.long <- try(cbind(Subject = 1:nrow(data), data) %>%
                       tidyr::gather(key = Day, value = Value, -Subject) %>%
                       arrange(Subject) %>%
                       mutate(Subject = as.factor(Subject),
                              Day = as.numeric({ gsub("Day", "", .$Day) })))
    
    fit.REML <- try(lmer(Value ~ 1 + Day + (1 + Day | Subject),
                         data = data.long, REML = TRUE), silent = TRUE)
    
    if(!inherits(fit.REML, "try-error")){
      
      out.REML <- est.REML(fit.REML)

      conv <- ifelse(any(grepl("failed to converge",
                               fit.REML@optinfo$conv$lme4$messages)), 0, 1)
      
      ci.norm.REML <- c("LB" = out.REML$est - 1.96*out.REML$se,
                        "UB" = out.REML$est + 1.96*out.REML$se)
      
      est[nsim*11 + i, ] <- c("jobid" = jobid, 
                              "iteration" = i,
                              "seed" = seed,
                              "method" = est.methods["REML"],
                              "ci.type" = conf.methods["REML"],
                              "tryerror" = 0,
                              "convergence" = conv,
                              "nobs" = nobs,
                              out.REML$est,
                              out.REML$se,
                              ci.norm.REML,
                              rep(NA, times = length(coef(fit))))  
      
    } else {
      
      est[nsim*11 + i, ] <- c("jobid" = jobid, 
                              "iteration" = i,
                              "seed" = seed,
                              "method" = est.methods["REML"],
                              "ci.type" = conf.methods["REML"],
                              "tryerror" = 1,
                              "convergence" = NA,
                              "nobs" = nobs,
                              rep(NA, times = length(all.info$estimates)-8))    
    }
    
    cat("Finished: i = ", i, "\t --- REML", "\n", sep = "")
    
  } else {
    
    cat("failure: i = ", i, "\t MLB: empty rows", "\n", sep = "")
    
    # Maximum Likelihood Bounds (MLB)
    est[nsim*0 + i, ] <- c("job.id" = jobid, 
                           "iteration" = i,
                           "seed" = seed,
                           "method" = est.methods["MLB"],
                           "ci.type" = conf.methods["normal"],
                           rep(NA, times = length(info)))
    
    # Jackknife Bounds: bias-corrected estimates (JB)
    a <- (i-1)*nobs + 1
    b <- i*nobs
    jack.resamples[a:b, ] <- cbind("jobid" = jobid, 
                                   "iteration" = i, 
                                   "seed" = seed,
                                   "resample" = 1:nobs,
                                   "tryerror" = 1,
                                   matrix(NA, 
                                          nrow = nobs, 
                                          ncol = length(all.info$resamples)-5))
    
    est[nsim*1 + i, ] <- c("jobid" = jobid, 
                           "iteration" = i,
                           "seed" = seed,
                           "method" = est.methods["JB"],
                           "ci.type" = conf.methods["JB"],
                           "tryerror" = 1,
                           "convergence" = NA,
                           "nobs" = nobs,
                           rep(NA, times = length(all.info$estimates)-8))  
    
    # Bootstrap Bounds: bias-corrected estimates (BB) - ci.norm
    a <- (i-1)*ndraws + 1
    b <- i*ndraws
    boot.resamples[a:b, ] <- cbind("jobid" = jobid, 
                                   "iteration" = i, 
                                   "seed" = seed,
                                   "resample" = 1:ndraws,
                                   "tryerror" = 1,
                                   matrix(NA, 
                                          nrow = ndraws, 
                                          ncol = length(all.info$resamples)-5))
    
    est[nsim*2 + i, ] <- c("jobid" = jobid, 
                           "iteration" = i,
                           "seed" = seed,
                           "method" = est.methods["BB"],
                           "ci.type" = conf.methods["BB-ci.norm"],
                           "tryerror" = 1,
                           "convergence" = NA,
                           "nobs" = nobs,
                           rep(NA, times = length(all.info$estimates)-8))   
    
    # Bootstrap Bounds: bias-corrected estimates (BB) - ci.basic
    est[nsim*3 + i, ] <- c("jobid" = jobid, 
                           "iteration" = i,
                           "seed" = seed,
                           "method" = est.methods["BB"],
                           "ci.type" = conf.methods["BB-ci.basic"],
                           "tryerror" = 1,
                           "convergence" = NA,
                           "nobs" = nobs,
                           rep(NA, times = length(all.info$estimates)-8)) 
    
    # Bootstrap Bounds: bias-corrected estimates (BB) - ci.stud
    est[nsim*4 + i, ] <- c("jobid" = jobid, 
                           "iteration" = i,
                           "seed" = seed,
                           "method" = est.methods["BB"],
                           "ci.type" = conf.methods["BB-ci.stud"],
                           "tryerror" = 1,
                           "convergence" = NA,
                           "nobs" = nobs,
                           rep(NA, times = length(all.info$estimates)-8))   
    
    # Bootstrap Bounds: bias-corrected estimates (BB) - ci.perc
    est[nsim*5 + i, ] <- c("jobid" = jobid, 
                           "iteration" = i,
                           "seed" = seed,
                           "method" = est.methods["BB"],
                           "ci.type" = conf.methods["BB-ci.perc"],
                           "tryerror" = 1,
                           "convergence" = NA,
                           "nobs" = nobs,
                           rep(NA, times = length(all.info$estimates)-8)) 
    
    # Bootstrap Bounds: bias-corrected estimates (BB) - ci.bca
    est[nsim*6 + i, ] <- c("jobid" = jobid, 
                           "iteration" = i,
                           "seed" = seed,
                           "method" = est.methods["BB"],
                           "ci.type" = conf.methods["BB-ci.bca"],
                           "tryerror" = 1,
                           "convergence" = NA,
                           "nobs" = nobs,
                           rep(NA, times = length(all.info$estimates)-8))  
    
    # Ozenne bias-corrected estimates (Ozenne)
    est[nsim*7 + i, ] <- c("jobid" = jobid, 
                           "iteration" = i,
                           "seed" = seed,
                           "method" = est.methods["Ozenne"],
                           "ci.type" = conf.methods["Ozenne"],
                           "tryerror" = 1,
                           "convergence" = NA,
                           "nobs" = nobs,
                           rep(NA, times = length(all.info$estimates)-8))
    
    # Simple bias-corrected estimates (simple): n/n-1
    est[nsim*8 + i, ] <- c("jobid" = jobid, 
                           "iteration" = i,
                           "seed" = seed,
                           "method" = est.methods["simple.min1"],
                           "ci.type" = conf.methods["simple.min1"],
                           "tryerror" = 1,
                           "convergence" = NA,
                           "nobs" = nobs,
                           rep(NA, times = length(all.info$estimates)-8))
    
    # Simple bias-corrected estimates (simple): n/n-2
    est[nsim*9 + i, ] <- c("jobid" = jobid, 
                           "iteration" = i,
                           "seed" = seed,
                           "method" = est.methods["simple.min2"],
                           "ci.type" = conf.methods["simple.min2"],
                           "tryerror" = 1,
                           "convergence" = NA,
                           "nobs" = nobs,
                           rep(NA, times = length(all.info$estimates)-8))
    
    # Simple bias-corrected estimates (simple): n/n-3
    est[nsim*10 + i, ] <- c("jobid" = jobid, 
                            "iteration" = i,
                            "seed" = seed,
                            "method" = est.methods["simple.min3"],
                            "ci.type" = conf.methods["simple.min3"],
                            "tryerror" = 1,
                            "convergence" = NA,
                            "nobs" = nobs,
                            rep(NA, times = length(all.info$estimates)-8))
    
    # REML
    est[nsim*11 + i, ] <- c("jobid" = jobid, 
                            "iteration" = i,
                            "seed" = seed,
                            "method" = est.methods["REML"],
                            "ci.type" = conf.methods["REML"],
                            "tryerror" = 1,
                            "convergence" = NA,
                            "nobs" = nobs,
                            rep(NA, times = length(all.info$estimates)-8))
    
  }
  
}

(finishing <- Sys.time())
cat("Finished ", nsim, " simulations! \n")
print(finishing - starting)
