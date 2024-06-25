# YR 25 June 2024
# small script to illustrate how to get the needed ingredients to compute
# the explicit bias reduction correction term

library(lavaan)
library(numDeriv)

# using the Holzinger & Swineford example, but any model would do
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
# fit the model
fit <- cfa(HS.model, data = HolzingerSwineford1939)
# ML estimates
coef(fit)

# starting values only (no optimization)
fit0 <- cfa(HS.model, data = HolzingerSwineford1939, do.fit = FALSE)
theta.init <- coef(fit0)

# extract some internal slots for housekeeping
lavsamplestats <- fit0@SampleStats
lavmodel       <- fit0@Model
lavdata        <- fit0@Data
lavoptions     <- fit0@Options

# extract ingredients needed for explicit bias reduction
# starting from the loglikelihood, not the objective function used
# by lavaan 

# theta = the parameter vector

# return loglikelihood (for a single group only), give a parameter vector 'theta'
obj_loglik <- function(theta) {
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

# numerical gradient (just to double-check)
grad.num <- numDeriv::grad(func = obj_loglik, x = theta.init)

# gradient
grad_loglik <- function(theta) {
  # fill in parameters in lavaan's internal matrix representation
  this.lavmodel <- lav_model_set_parameters(lavmodel, x = theta)
  # gradient of F_ML (not loglik yet) 
  grad.F <- lavaan:::lav_model_gradient(lavmodel = this.lavmodel, 
              lavsamplestats = lavsamplestats, lavdata = lavdata) 
  # rescale so we get gradient of loglik
  N <- lavsamplestats@ntotal
  grad.loglik <- -1 * N * grad.F
  grad.loglik
}

# at this point, we can use quasi-newton to find the ML estimates of theta;
# but note that nlminb() only minimizes, so we need to change the sign
min_obj_loglik <- function(theta) {
  loglik <- obj_loglik(theta)
  # minimize
  -1 * loglik
}
min_grad_loglik <- function(theta) {
  grad_loglik <- grad_loglik(theta)
  # minimize
  -1 * grad_loglik
}
est1 <- nlminb(start = theta.init, objective = min_obj_loglik, 
               gradient = min_grad_loglik)
round(est1$par, 3)



# nummeric hessian (just to double-check)
hessian.num1 <- numDeriv::hessian(func = obj_loglik, x = theta.init)
# better approach: just Jacobian of the (analytic) gradient
hessian.num <- numDeriv::jacobian(func = grad_loglik, x = theta.init)

# hessian
hessian_loglik <- function(theta) {
  # fill in parameters in lavaan's internal matrix representation
  this.lavmodel <- lav_model_set_parameters(lavmodel, x = theta)
  # gradient of F_ML (not loglik yet) 
  hessian.F <- lavaan:::lav_model_hessian(lavmodel = this.lavmodel,
              lavsamplestats = lavsamplestats, lavdata = lavdata,
              lavoptions = lavoptions)
  # rescale so we get gradient of loglik
  N <- lavsamplestats@ntotal
  hessian.loglik <- -1 * N * hessian.F
  hessian.loglik
}

# at this point, we can use newton-raphson to find the ML estimates of theta;
# but note that nlminb() only minimizes, so we need to change the sign
min_hessian_loglik <- function(theta) {
  hessian.loglik <- hessian_loglik(theta)
  -1 * hessian.loglik
}
est2 <- nlminb(start = theta.init, objective = min_obj_loglik, 
               gradient = min_grad_loglik, hessian = min_hessian_loglik)
# only 9 iterations were needed, but each iteration is (much) slower
# as it takes time to compute the hessian every time
round(est2$par, 3)


# casewise scores (score = grad of loglik for a single observation)
scores_loglik <- function(theta) {
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

# outer product of casewise scores, divided by N
first_order_unit_information_loglik <- function(theta) {
  # ill in parameters in lavaan's internal matrix representation
  this.lavmodel <- lav_model_set_parameters(lavmodel, x = theta)
  # compute model implied mean and (co)variance matrix
  lavimplied <- lav_model_implied(this.lavmodel)

  info <- lavaan:::lav_model_information_firstorder(
    lavmodel = this.lavmodel, lavdata = lavdata, lavsamplestats = lavsamplestats,
    lavoptions = lavoptions
  )
  # no scaling needed
  info
}

# compare
# tmp1 <- crossprod(scores_loglik(theta.init))/nobs(fit0)
# versus
# tmp2 <- first_order_unit_information_loglik(theta.init)
# range(tmp1 - tmp2)

# ingredients for bias reduction
j_theta <- function(theta) {
  j.theta <- hessian_loglik(theta) / nobs(fit0) # I think? same scale as e_theta?
  j.theta
}

# variability matrix
e_theta <- function(theta) {
  e.theta <- first_order_unit_information_loglik(theta)
  e.theta
}

A_theta <- function(theta) {
  ff <- function(theta) {
    j.theta <- j_theta(theta)
    j.theta.inv <- solve(j.theta)
    e.theta <- e_theta(theta)
    #trace_theta <- sum(diag(j.theta.inv %*% e.theta))
    trace_theta <- sum(j.theta.inv * e.theta) # same as trace
    -0.5 * trace_theta
  }
  out <- numDeriv::grad(func = ff, x = theta)
  out
}

# ML estimates
theta.ml <- est2$par

# bias reduction
j.theta.ml <- j_theta(theta.ml)
A.theta.ml <- A_theta(theta.ml)

theta.correction <- j.theta.ml %*% A.theta.ml
# hm? too big?
# does this need rescaling?
# divide by N perhaps?

# once ok, then:
# theta.ml.br <- theta.ml + theta.correction 

