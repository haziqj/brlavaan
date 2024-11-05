##==============================================================================
## Project:     Resampling based bias correction for small sample SEM
##  
## Script:      Simulations two-factor model
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
                 "simple.min3" = 7)

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
                  "simple.min3" = 11)
                  
#-------------------------------------------------------------------------------
# Load libraries
#-------------------------------------------------------------------------------

library(lavaan)
library(dplyr)

#-------------------------------------------------------------------------------
# Source files
#-------------------------------------------------------------------------------

source("00_Functions/00_Function_sd_lav.extract.vector.R")
source("00_Functions/00_Function_sd_jack_draws.R")
source("00_Functions/00_Function_sd_boot_draws.R")
source("00_Functions/00_Function_sd_jack_bc.R")
source("00_Functions/00_Function_sd_boot_bc.R")
source("00_Functions/00_Function_oz_bc_2FSEM.R")
source("00_Functions/00_Function_sd_varcomp_rescaling.R")

#-------------------------------------------------------------------------------
# Set up data generating model
#-------------------------------------------------------------------------------

n.factors <- 2
n.indicators <- 6
  
# Create (empty) matrices
lambda <- matrix(data = 0, nrow = n.indicators, ncol = n.factors)
psi <- matrix(data = 0, nrow = n.factors, ncol = n.factors)
beta <- matrix(data = 0, nrow = n.factors, ncol = n.factors)
theta <- matrix(data = 0, nrow = n.indicators, ncol = n.indicators)
eta <- matrix(data = NA, nrow = n.factors, ncol = nobs)
e <- matrix(data = NA, nrow = n.indicators, ncol = nobs)
  
# Define aspired reliability of indicators
reliability <- matrix(data = 0, nrow = n.indicators, ncol = n.indicators)
diag(reliability) <- readr::parse_number(rel)/100
  
# Define factor loadings, factor variances and residual variances
lambda[1:6, 1] <- c(1, 0.7, 0.6, 0, 0, 0)
lambda[1:6, 2] <- c(0, 0, 0, 1, 0.7, 0.6)
diag(psi) <- c(1, 1)
beta[2, 1] <- 0.25
diag(theta) <- diag((lambda %*% psi %*% t(lambda)) %*% 
                      solve(reliability) - (lambda %*% psi %*% t(lambda)))

# Population model in lavaan syntax
if(readr::parse_number(rel)/100 == 0.8) {
  
  pop.model <- 'fx =~ 1*x1 + 0.7*x2 + 0.6*x3
                fy =~ 1*y1 + 0.7*y2 + 0.6*y3
                
                fy ~ 0.25*fx
                
                x1 ~~ 0.25*x1
                x2 ~~ 0.1225*x2
                x3 ~~ 0.09*x3
                
                y1 ~~ 0.25*y1
                y2 ~~ 0.1225*y2
                y3 ~~ 0.09*y3
               '
  
} else if(readr::parse_number(rel)/100 == 0.5) {
  
  pop.model <- 'fx =~ 1*x1 + 0.7*x2 + 0.6*x3
                fy =~ 1*y1 + 0.7*y2 + 0.6*y3
                
                fy ~ 0.25*fx
                
                x1 ~~ 1*x1
                x2 ~~ 0.49*x2
                x3 ~~ 0.36*x3
                
                y1 ~~ 1*y1
                y2 ~~ 0.49*y2
                y3 ~~ 0.36*y3
               '
} else {cat("'rel' needs to be 0.5 or 0.8")}

# Specify analysis model
model <- 'fx =~ x1 + x2 + x3
          fy =~ y1 + y2 + y3
          fy ~ fx'

#-------------------------------------------------------------------------------
# Prepare containers to store estimates
#-------------------------------------------------------------------------------

# Test model fit (shortcut to extract all info)
fit <- cfa(data = simulateData(pop.model, empirical = TRUE, sample.nobs = 15), 
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
    
    random <- covsim::rIG(N = nobs, 
                          sigma.target = psi, 
                          skewness = rep(0, times = ncol(lambda)),
                          excesskurtosis = rep(0, times = ncol(lambda)),
                          reps = 1)[[1]]
    
    fx <- random[ , 1]
    fy <- fx*beta[2, 1] + random[ , 2]
    factors <- rbind(fx, fy)
    
    e <- t(covsim::rIG(N = nobs, 
                       sigma.target = theta, 
                       skewness = rep(0, times = nrow(theta)), 
                       excesskurtosis = rep(0, times = nrow(theta)), 
                       reps = 1)[[1]])
    
  } else if(dist == "Kurtosis") {
    
    random <- covsim::rIG(N = nobs, 
                          sigma.target = psi, 
                          skewness = rep(0, times = ncol(lambda)),
                          excesskurtosis = rep(6, times = ncol(lambda)),
                          reps = 1)[[1]]
    
    fx <- random[ , 1]
    fy <- fx*beta[2, 1] + random[ , 2]
    factors <- rbind(fx, fy)
    
    e <- t(covsim::rIG(N = nobs, 
                       sigma.target = theta, 
                       skewness = rep(0, times = nrow(theta)), 
                       excesskurtosis = rep(6, times = nrow(theta)), 
                       reps = 1)[[1]])
    
  } else if(dist == "NonNormal") {
    
    random <- covsim::rIG(N = nobs, 
                          sigma.target = psi, 
                          skewness = rep(-2, times = ncol(lambda)),
                          excesskurtosis = rep(6, times = ncol(lambda)),
                          reps = 1)[[1]]
    
    fx <- random[ , 1]
    fy <- fx*beta[2, 1] + random[ , 2]
    factors <- rbind(fx, fy)
    
    e <- t(covsim::rIG(N = nobs, 
                       sigma.target = theta, 
                       skewness = rep(-2, times = nrow(theta)), 
                       excesskurtosis = rep(6, times = nrow(theta)), 
                       reps = 1)[[1]])
    
  } else {cat("'dist' needs to be 'Normal', 'Kurtosis' or 'NonNormal'")}

  data <- as.data.frame(t(lambda %*% factors + e))
  colnames(data) <- c("x1", "x2", "x3", "y1", "y2", "y3")
  
  #-----------------------------------------------------------------------------
  # Fit model
  #-----------------------------------------------------------------------------
  
  fit <- try(cfa(data = data, 
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
    ozenne.out <- try(yr_ozenne_bcs(fit, return.se = TRUE), silent = TRUE)
    
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
    
  }
    
}

(finishing <- Sys.time())
cat("Finished ", nsim, " simulations! \n")
print(finishing - starting)
