# Log-likelihood function
loglik <- function(
    theta,  # packed / short version
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
  loglik
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
  N <- lavsamplestats@ntotal
  # hessian.loglik <- -1 * N * hessian.F
  H <- N * hessian.F

  # repack hessian
  if (lavmodel@eq.constraints) {
    H <- t(lavmodel@eq.constraints.K) %*% H %*% lavmodel@eq.constraints.K
  }

  H
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

  # repack
  if (lavmodel@eq.constraints) {
    score_matrix <- score_matrix %*% lavmodel@eq.constraints.K
  }

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
  out <- info * lavsamplestats@ntotal

  # repack
  if (lavmodel@eq.constraints) {
    out <- t(lavmodel@eq.constraints.K) %*% out %*% lavmodel@eq.constraints.K
  }

  out
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
  j <- hessian_loglik(
    theta = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
  jinv <- try(solve(j), silent = TRUE)

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
  jinv <- try(solve(j), silent = TRUE)
  A <- numDeriv::grad(
    func = penalty,
    x = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
  -drop(jinv %*% A)
}

fit_sem <- function(
    model,
    data,
    method = "ML",
    debug = FALSE,
    theta_init = NULL
  ) {

  method   <- match.arg(method, c("ML", "iRBM", "eRBM", "iRBMp"))
  is_ML    <- method == "ML"
  is_iRBM  <- method == "iRBM"
  is_eRBM  <- method == "eRBM"
  is_iRBMp <- method == "iRBMp"

  # Initialise {lavaan} model object
  fit0 <- sem(model = model, data = data, do.fit = FALSE)
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
        loglik = loglik(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions),
        grad_loglik = grad_loglik(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions),
        j = hessian_loglik(theta_pack, lavmodel, lavsamplestats, lavdata, lavoptions),
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
      lavoptions = lavoptions
    )
    b <- 0
    if (is_eRBM) b <- b + bias(
      theta = res$par,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
    est <- res$par - b
  } else {
    if (is_iRBMp) {
      # Implicit RBM with plugin penalty
      obj <- function(x, ...) -loglik(x, ...) - penalty(x, ...) + sum(x ^ 2) / n
    } else {
      # Implicit RBM
      obj <- function(x, ...) -loglik(x, ...) - penalty(x, ...)
    }
    res <- nlminb(
      start = theta_pack,
      objective = obj,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions
    )
    est <- res$par
  }
  end_time <- Sys.time()

  # Standard errors ------------------------------------------------------------
  j <- hessian_loglik(
    theta = as.numeric(est),
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions
  )
  vars <- try(diag(solve(j)), silent = TRUE)
  if (inherits(vars, "try-error")) vars <- rep(NA, length(est))

  # unpack est and se
  if (lavmodel@eq.constraints) {
    est  <- as.numeric(lavmodel@eq.constraints.K %*% est) +
      lavmodel@eq.constraints.k0
    vars <- as.numeric(lavmodel@eq.constraints.K %*% vars)
  }

  list(
    coefficients = est,
    stderr = sqrt(vars),
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
