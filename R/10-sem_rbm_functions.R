# Log-likelihood function
loglik <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions,
    idx
  ) {

  theta <- expand_theta(theta, idx)
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
    lavoptions,
    idx
  ) {
  theta <- expand_theta(theta, idx)
  # fill in parameters in lavaan's internal matrix representation
  this.lavmodel <- lav_model_set_parameters(lavmodel, x = theta)
  # gradient of F_ML (not loglik yet)
  grad.F <- lavaan:::lav_model_gradient(lavmodel = this.lavmodel,
                                        lavsamplestats = lavsamplestats, lavdata = lavdata)
  # rescale so we get gradient of loglik
  N <- lavsamplestats@ntotal
  grad.loglik <- -1 * N * grad.F
  grad.loglik[unique(idx)]
}

# Hessian
hessian_loglik <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions,
    idx
  ) {
  theta <- expand_theta(theta, idx)
  # fill in parameters in lavaan's internal matrix representation
  this.lavmodel <- lav_model_set_parameters(lavmodel, x = theta)
  # gradient of F_ML (not loglik yet)
  hessian.F <- lavaan:::lav_model_hessian(lavmodel = this.lavmodel,
                                          lavsamplestats = lavsamplestats, lavdata = lavdata,
                                          lavoptions = lavoptions)
  # rescale so we get gradient of loglik
  N <- lavsamplestats@ntotal
  # hessian.loglik <- -1 * N * hessian.F
  hessian.loglik <- N * hessian.F  # FIXME: Multiply by -1 or not?
  hessian.loglik[unique(idx), unique(idx)]
}

# Casewise scores (score = grad of loglik for a single observation)
scores_loglik <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions,
    idx
  ) {
  theta <- expand_theta(theta, idx)
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
  score_matrix[, unique(idx)]
}

# Outer product of casewise scores
first_order_unit_information_loglik <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions,
    idx
  ) {
  theta <- expand_theta(theta, idx)
  # fill in parameters in lavaan's internal matrix representation
  this.lavmodel <- lav_model_set_parameters(lavmodel, x = theta)
  # compute model implied mean and (co)variance matrix
  lavimplied <- lav_model_implied(this.lavmodel)

  info <- lavaan:::lav_model_information_firstorder(
    lavmodel = this.lavmodel, lavdata = lavdata, lavsamplestats = lavsamplestats,
    lavoptions = lavoptions
  )

  # this is unit information, so multiply by N
  out <- info * lavsamplestats@ntotal
  out[unique(idx), unique(idx)]
}

penalty <- function(
    theta,
    lavmodel,
    lavsamplestats,
    lavdata,
    lavoptions,
    idx
  ) {
  e <- first_order_unit_information_loglik(
    theta = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    idx = idx
  )
  j <- hessian_loglik(
    theta = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    idx = idx
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
    lavoptions,
    idx
  ) {
  j <- hessian_loglik(
    theta = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    idx = idx
  )
  jinv <- try(solve(j), silent = TRUE)
  A <- numDeriv::grad(
    func = penalty,
    x = theta,
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    idx = idx
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

  # Parameter constraints
  idx_constr <- get_constr_idx(pt)
  theta_init <- contract_theta(theta_init, idx_constr)

  # DEBUG ----------------------------------------------------------------------
  if (isTRUE(debug)) {
    out <- with(
      list(
        theta = theta_init,
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats,
        lavdata = lavdata,
        lavoptions = lavoptions,
        idx = idx_constr
      ),
      list(
        theta_nlimb = theta_init,
        theta = expand_theta(theta_init, idx),
        n = n,
        loglik = loglik(theta, lavmodel, lavsamplestats, lavdata, lavoptions, idx),
        grad_loglik = grad_loglik(theta, lavmodel, lavsamplestats, lavdata, lavoptions, idx),
        j = hessian_loglik(theta, lavmodel, lavsamplestats, lavdata, lavoptions, idx),
        e = first_order_unit_information_loglik(theta, lavmodel, lavsamplestats, lavdata, lavoptions, idx),
        scores_loglik = scores_loglik(theta, lavmodel, lavsamplestats, lavdata, lavoptions, idx),
        penalty = penalty(theta, lavmodel, lavsamplestats, lavdata, lavoptions, idx),
        bias = bias(theta, lavmodel, lavsamplestats, lavdata, lavoptions, idx)
      )
    )
    return(out)
  }

  # Start fitting process ------------------------------------------------------
  start_time <- Sys.time()
  if (is_ML | is_eRBM) {
    # ML or explicit RBM -- either way requires MLE
    res <- nlminb(
      start = theta_init,
      objective = function(x, ...) -loglik(x, ...),
      gradient = function(x, ...) -grad_loglik(x, ...),
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      idx = idx_constr
    )
    b <- 0
    if (is_eRBM) b <- b + bias(
      theta = res$par,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      idx = idx_constr
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
      start = theta_init,
      objective = obj,
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      idx = idx_constr
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
    lavoptions = lavoptions,
    idx = idx_constr
  )
  sds <- try(sqrt(diag(solve(j))), silent = TRUE)
  if (inherits(sds, "try-error")) sds <- rep(NA, length(est))

  list(
    coefficients = expand_theta(est, idx_constr),
    stderr = expand_theta(sds, idx_constr),
    time = end_time - start_time
  )
}

get_lav_stuff <- function(fit) {
  list(
    lavmodel       = fit@Model,
    lavsamplestats = fit@SampleStats,
    lavdata        = fit@Data,
    lavoptions     = fit@Options,
    idx_constr     = get_constr_idx(partable(fit))
  )
}

get_constr_idx <- function(pt) {
  free_params <- plabs <- pt$plabel[pt$free > 0]
  idx <- seq_along(free_params)
  equal_pt <- pt[pt$op == "==", ]
  if (nrow(equal_pt) > 0) {
    for (i in seq_len(nrow(equal_pt))) {
      cur_rhs <- equal_pt$rhs[i]
      cur_lhs <- equal_pt$lhs[i]
      free_params[plabs == cur_rhs] <- cur_lhs
    }

    idx <- match(free_params, unique(free_params))
  }
  idx
}

contract_theta <- function(theta, idx) {
  if (is.null(idx))
    return(theta)
  else
    return(theta[unique(idx)])
}

expand_theta <- function(theta, idx) {
  if (is.null(idx))
    return(theta)
  else {
    return(theta[idx])
  }
}

