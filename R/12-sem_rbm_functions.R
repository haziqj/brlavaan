# copy of 10-sem_rbm_functions.R with the following changes:
# not for simulation purposes, only for debugging
#
# 1) loglik() has RBM and plugin.penalty arguments
# 2) we check explicitly for non-positive definite Sigma
# 3) we check for non-convergence (and return NAs)
# 4) we print out the loglik/rbm/plugin.penalty terms at each iteration
# 5) for iRBM and iRBMp, we use EXPECTED information (only for numerical reasons)

# Log-likelihood function
loglik <- function(
    theta,  # packed / short version
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions,
    RBM = FALSE,
    plugin.penalty = FALSE
  ) {

  # unpack theta
  theta.packed <- theta
  if (lavmodel@eq.constraints) {
    theta <- as.numeric(lavmodel@eq.constraints.K %*% theta) +
      lavmodel@eq.constraints.k0
  }

  # fill in parameters in lavaan's internal matrix representation
  this.lavmodel <- lav_model_set_parameters(lavmodel, x = theta)

  # compute model implied mean and (co)variance matrix
  lavimplied <- lav_model_implied(this.lavmodel)

  # check if lavimplied$cov[[1]] is PD
  eigvals <- eigen(lavimplied$cov[[1]], symmetric = TRUE,
                   only.values = TRUE)$values
  # cat("eigvals = ", eigvals, "\n")
  if (any(eigvals < 1e-07)) {
    # return huge number, to signal the nlminb() optimizer something is not
    # quite ok with this parameter vector
    return(-1e40)
  }

  # meanstructure?
  if (lavmodel@meanstructure) {
    Mu <- lavimplied$mean[[1]]
  } else {
    Mu <- lavsamplestats@mean[[1]]
  }

  # total loglik
  loglik <- lavaan:::lav_mvnorm_loglik_samplestats(
    sample.mean = lavsamplestats@mean[[1]],
    sample.cov  = lavsamplestats@cov[[1]],
    sample.nobs = lavsamplestats@nobs[[1]],
    Mu          = Mu,
    Sigma       = lavimplied$cov[[1]],
    x.idx       = lavsamplestats@x.idx[[1]],
    x.mean      = lavsamplestats@mean.x[[1]],
    x.cov       = lavsamplestats@cov.x[[1]],
    Sinv.method = "eigen",
    Sigma.inv   = NULL
  )

  # RBM
  rbm.term <- penalty.term <- 0
  if (RBM) {
    rbm.term <- penalty(theta = theta.packed,
                        lavmodel = lavmodel,
                        lavsamplestats = lavsamplestats,
                        lavdata = lavdata,
                        lavoptions = lavoptions)
  }
  if (plugin.penalty) {
    penalty.term <- sum(theta.packed^2) / lavsamplestats@ntotal
  }

  return.value <- loglik + rbm.term - penalty.term

  # nlminb() works better with smaller numbers
  #return.value <- return.value / (2 * lavsamplestats@ntotal)

  cat("loglik =", loglik, "rbm.term = ", rbm.term, "plugin pen = ",
      penalty.term, "\n")

  return.value
}

# Gradient function
grad_loglik <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions
  ) {

  # unpack theta
  if (lavmodel@eq.constraints) {
    theta <- as.numeric(lavmodel@eq.constraints.K %*% theta) +
      lavmodel@eq.constraints.k0
  }

  # fill in parameters in lavaan's internal matrix representation
  this.lavmodel <- lav_model_set_parameters(lavmodel, x = theta)
  # gradient of F_ML (not loglik yet)
  grad.F <- lavaan:::lav_model_gradient(lavmodel = this.lavmodel,
                                        lavsamplestats = lavsamplestats,
                                        lavdata = lavdata)
  # rescale so we get gradient of loglik
  N <- lavsamplestats@ntotal
  grad.loglik <- -1 * N * grad.F

  # repack gradients
  if (lavmodel@eq.constraints) {
    grad.loglik <- as.numeric(grad.loglik %*% lavmodel@eq.constraints.K)
  }

  # adjust gradient for rescaling
  # grad.loglik <- grad.loglik / (2 * lavsamplestats@ntotal)
  
  grad.loglik
}

# Hessian
hessian_loglik <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions
  ) {

  # unpack theta
  if (lavmodel@eq.constraints) {
    theta <- as.numeric(lavmodel@eq.constraints.K %*% theta) +
      lavmodel@eq.constraints.k0
  }

  # fill in parameters in lavaan's internal matrix representation
  this.lavmodel <- lav_model_set_parameters(lavmodel, x = theta)
  # gradient of F_ML (not loglik yet)
  hessian.F <- lavaan:::lav_model_hessian(lavmodel = this.lavmodel,
                                          lavsamplestats = lavsamplestats,
                                          lavdata = lavdata,
                                          lavoptions = lavoptions)

  # rescale so we get gradient of loglik
  lavsamplestats@ntotal * hessian.F  # FULL INFORMATION
}

# observed information 
observed_information_loglik <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions,
    simple = FALSE
  ) {

  # change to 'simplified' hessian
  if (simple) {
    lavoptions$observed.information <- c("h1", "h1")
  }
  
  # unpack theta
  if (lavmodel@eq.constraints) {
    theta <- as.numeric(lavmodel@eq.constraints.K %*% theta) +
      lavmodel@eq.constraints.k0
  } 
    
  # fill in parameters in lavaan's internal matrix representation
  this.lavmodel <- lav_model_set_parameters(lavmodel, x = theta)

  obs.info <- lavaan:::lav_model_information_observed(lavmodel = this.lavmodel,
                                          lavsamplestats = lavsamplestats,
                                          lavdata = lavdata,
                                          lavoptions = lavoptions)
  # this is unit information, so multiply by N
  lavsamplestats@ntotal * obs.info
}

# expected information 
expected_information_loglik <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions
  ) {
  
  # unpack theta
  if (lavmodel@eq.constraints) {
    theta <- as.numeric(lavmodel@eq.constraints.K %*% theta) +
      lavmodel@eq.constraints.k0
  } 
    
  # fill in parameters in lavaan's internal matrix representation
  this.lavmodel <- lav_model_set_parameters(lavmodel, x = theta)
    
  obs.info <- lavaan:::lav_model_information_expected(lavmodel = this.lavmodel,
                                          lavsamplestats = lavsamplestats,
                                          lavdata = lavdata,
                                          lavoptions = lavoptions)
  # this is unit information, so multiply by N
  lavsamplestats@ntotal * obs.info
} 




# Casewise scores (score = grad of loglik for a single observation)
scores_loglik <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions
  ) {

  # unpack theta
  if (lavmodel@eq.constraints) {
    theta <- as.numeric(lavmodel@eq.constraints.K %*% theta) +
      lavmodel@eq.constraints.k0
  }

  # fill in parameters in lavaan's internal matrix representation
  this.lavmodel <- lav_model_set_parameters(lavmodel, x = theta)
  # compute model implied mean and (co)variance matrix
  lavimplied <- lav_model_implied(this.lavmodel)
  # only 1 group
  moments <- list(cov = lavimplied$cov[[1]])
  if(lavmodel@meanstructure) {
    moments$mean <- lavimplied$mean[[1]]
  }

  ntab <- unlist(lavdata@norig)
  ntot <- sum(ntab)
  npar <- length(theta)
  score_matrix <- lavaan:::lav_scores_ml(
    ntab = ntab, ntot = ntot, npar = npar,
    moments = moments, lavdata = lavdata, lavsamplestats = lavsamplestats,
    lavmodel = lavmodel, lavoptions = lavoptions, scaling = FALSE
  )

  score_matrix
}

# Outer product of casewise scores
first_order_unit_information_loglik <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions
  ) {

  # unpack theta
  if (lavmodel@eq.constraints) {
    theta <- as.numeric(lavmodel@eq.constraints.K %*% theta) +
      lavmodel@eq.constraints.k0
  }

  # fill in parameters in lavaan's internal matrix representation
  this.lavmodel <- lav_model_set_parameters(lavmodel, x = theta)
  # compute model implied mean and (co)variance matrix
  lavimplied <- lav_model_implied(this.lavmodel)

  info <- lavaan:::lav_model_information_firstorder(
    lavmodel = this.lavmodel,
    lavdata = lavdata,
    lavsamplestats = lavsamplestats,
    lavoptions = lavoptions
  )

  # this is unit information, so multiply by N
  lavsamplestats@ntotal * info
}

penalty <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions
  ) {
  e <- first_order_unit_information_loglik(
    theta = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
  j <- expected_information_loglik(
    theta = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
  if (lavmodel@eq.constraints) {
    jinv <- lavaan:::lav_model_information_augment_invert(
      lavmodel = lavmodel,
      information = j,
      inverted = TRUE,
      check.pd = FALSE,
      use.ginv = FALSE,
      rm.idx = integer(0L)
    )
  } else {
    jinv <- solve(j)
  }

  if (inherits(jinv, "try-error")) {
    return(NA)
  } else {
    return(-sum(jinv * e) / 2)
  }
}

bias <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions
  ) {
  j <- hessian_loglik(
    theta = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
  jinv <- lavaan:::lav_model_information_augment_invert(
    lavmodel = lavmodel,
    information = j,
    inverted = TRUE,
    check.pd = FALSE,
    use.ginv = FALSE,
    rm.idx = integer(0L)
  )
  A <- numDeriv::grad(
    func = penalty,
    x = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )

  # unpack A
  if (lavmodel@eq.constraints) {
    A <- lavmodel@eq.constraints.K %*% A
  }

  out <- -drop(jinv %*% A)

  # repack! (since optimiser expects packed version)
  if (lavmodel@eq.constraints) {
    out <- as.numeric(out %*% lavmodel@eq.constraints.K)
  }
  out
}

fit_sem <- function(
    model,
    data,
    method = "ML",
    debug = FALSE,
    theta_init = NULL,
    trace = 0,
    lavfun = "sem",
    meanstructure = FALSE,
    bounds = "none"
  ) {

  method   <- match.arg(method, c("ML", "iRBM", "eRBM", "iRBMp"))
  is_ML    <- method == "ML"
  is_iRBM  <- method == "iRBM"
  is_eRBM  <- method == "eRBM"
  is_iRBMp <- method == "iRBMp"

  # Initialise {lavaan} model object
  lavargs <- list(
    model = model,
    data = data,
    do.fit = FALSE,
    meanstructure = meanstructure,
    bounds = bounds
  )
  fit0           <- do.call(lavfun, lavargs)
  lavmodel       <- fit0@Model
  lavsamplestats <- fit0@SampleStats
  lavdata        <- fit0@Data
  lavoptions     <- fit0@Options
  pt             <- partable(fit0)

  if (is.null(theta_init)) theta_init <- coef(fit0)  # starting values
  n <- nobs(fit0)

  # Pack theta
  if (lavmodel@eq.constraints) {
    theta_pack <- as.numeric(
      (theta_init - lavmodel@eq.constraints.k0) %*% lavmodel@eq.constraints.K
    )
  } else {
    theta_pack <- theta_init
  }

  # DEBUG ----------------------------------------------------------------------
  if (isTRUE(debug)) {
    out <- with(
      list(
        theta = theta_init,
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats,
        lavdata = lavdata,
        lavoptions = lavoptions
      ),
      list(
        theta_unpack = theta_init,
        theta_pack = theta_pack,
        n = n,
        lavmodel_eq.constraints = lavmodel@eq.constraints,
        lavmodel_eq.constraints.k0 = lavmodel@eq.constraints.k0,
        lavmodel_eq.constraints.K = lavmodel@eq.constraints.K,
        loglik = loglik(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions),
        grad_loglik = grad_loglik(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions),
        j = hessian_loglik(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions),
        jinv = lavaan:::lav_model_information_augment_invert(
          lavmodel = lavmodel,
          information = hessian_loglik(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions),
          inverted = TRUE,
          check.pd = FALSE,
          use.ginv = FALSE,
          rm.idx = integer(0L)
        ),
        e = first_order_unit_information_loglik(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions),
        scores_loglik = scores_loglik(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions),
        penalty = penalty(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions),
        bias = bias(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions)
      )
    )
    return(out)
  }

  # Start fitting process ------------------------------------------------------
  start_time <- Sys.time()
  if (is_ML | is_eRBM) {
    # ML or explicit RBM -- either way requires MLE
    res <- nlminb(
      start = theta_pack,
      objective = function(x, ...) -loglik(x, ...),
      gradient = function(x, ...) -grad_loglik(x, ...),
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      control = list(trace = trace)
    )
    # check convergence
    converged.flag <- TRUE
    b <- 0
    if (res$convergence != 0L) {
      converged.flag <- FALSE
      est <- rep(as.numeric(NA), length(res$par))
    } else {
      if (is_eRBM) b <- b + bias(
        theta = res$par,
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats,
        lavdata = lavdata,
        lavoptions = lavoptions
      )
      est <- res$par - b
    }
  } else {
    if (is_iRBMp) {
      # Implicit RBM with plugin penalty
      obj <- function(x, ...) -loglik(x, ..., RBM = TRUE, plugin.penalty = TRUE)
    } else {
      # Implicit RBM
      obj <- function(x, ...) -loglik(x, ..., RBM = TRUE)
    }
    res <- nlminb(
      start = theta_pack,
      objective = obj,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      control = list(trace = trace)
    )
    if (res$convergence != 0L) {
      converged.flag <- FALSE
      est <- rep(as.numeric(NA), length(res$par))
    } else {
      est <- res$par
    }
  }
  end_time <- Sys.time()

  # Standard errors ------------------------------------------------------------
  # j <- hessian_loglik(
  #   theta = as.numeric(est),
  #   lavmodel = lavmodel,
  #   lavsamplestats = lavsamplestats,
  #   lavdata = lavdata,
  #   lavoptions = lavoptions
  # )
  # jinv <- lavaan:::lav_model_information_augment_invert(
  #   lavmodel = lavmodel,
  #   information = j,
  #   inverted = TRUE,
  #   check.pd = FALSE,
  #   use.ginv = FALSE,
  #   rm.idx = integer(0L)
  # )
  # sds <- sqrt(diag(jinv))
  sds <- rep(as.numeric(NA), length(est))

  # unpack est
  if (lavmodel@eq.constraints) {
    est  <- as.numeric(lavmodel@eq.constraints.K %*% est) +
      lavmodel@eq.constraints.k0
  }

  list(
    coefficients = est,
    stderr = sds,
    time = end_time - start_time
  )
}

get_lav_stuff <- function(fit) {
  # Utility function to extract lavaan stuff
  list(
    lavmodel       = fit@Model,
    lavsamplestats = fit@SampleStats,
    lavdata        = fit@Data,
    lavoptions     = fit@Options
  )
}
