# This file illustrates two ways to compute the Hessian matrix
# for the Holzinger & Swineford example (a 3-factor CFA)
# 1) numerically (using numDeriv)
# 2) analytically (using eq 15 of Savalei & Rosseel, 2022)

library(lavaan)
library(numDeriv)

Data <- HolzingerSwineford1939[, paste0("x",1:9)] 
N <- nrow(Data)

#############
#### H1 #####
#############

# saturated model (h1)
h1.Mu <- colMeans(Data)
h1.Sigma <- cov(Data) * (N-1)/N # ML estimate


#############
#### H0 #####
#############

model <- '
    visual  =~ x1 + x2 + x3
    textual =~ x4 + x5 + x6
    speed   =~ x7 + x8 + x9

    x1 ~ start(1)*1
    x2 ~ start(1)*1
    x3 ~ start(1)*1
    x4 ~ start(1)*1
    x5 ~ start(1)*1
    x6 ~ start(1)*1
    x7 ~ start(1)*1
    x8 ~ start(1)*1
    x9 ~ start(1)*1
'
# we stop the iterations too early, so the solution is NOT yet the ML solution
# this is just to get a reasonable (but not ML) model implied estimates
# for Mu and Sigma
fit <- sem(model, data = Data, estimator = "ML", meanstructure = TRUE,
           optim.attempts = 1L,
           control = list(iter.max = 10L), verbose = TRUE)
implied <- lavTech(fit, "implied")[[1]]
h0.Mu    <- implied$mean
h0.Sigma <- implied$cov

# not the same, so meanstrucutre is not saturated
h0.Mu - h1.Mu

# this is what lavaan does by default:
# using a numerical approximation of the jacobian of the analytic gradient:
A.theta.num <- lavTech(fit, "Hessian")
# note: this is not the Hessian of the logl, but of the discrepancy function
# F.ML!!
# as a result, it is scaled by -1/N 
# (it really should be -2/N, but currently lavaan still divides F by 2)

# order of Hessian corresponds to coef:
# - 6 factor loadings
# - 9 (residual) ov variances 
# - 3 lv variances
# - 3 lv covariances
# - 9 ov intercepts
coef(fit)

# compute Delta = d c(Mu, vech(Sigma)) / d theta
# - columns == 30 free parameters
# - 1-9 rows Mu
# - 10-54 rows: vech(Sigma)
Delta <- lavTech(fit, "Delta")[[1]]

# using h1.hessian.analytic, but use the h0 Mu and Sigma
# - first 9 rows/cols: Mu
# - rest: vech(Sigma)
A.beta <-
    -1/N * lavaan:::lav_mvnorm_logl_hessian_samplestats(sample.mean = h1.Mu, 
        sample.cov = h1.Sigma, sample.nobs = N, Mu = h0.Mu, Sigma = h0.Sigma)

# using 'simplified' Hessian -- see eq (16) in Savalei & Rosseel (2022)
A.theta.h1 <- t(Delta) %*% A.beta %*% Delta

# this is not the same as the 'real' Hessian, but close:
range(A.theta.num - A.theta.h1)
# and most entries are (almost) zero!
zapsmall(A.theta.num - A.theta.h1)


# compute the 'H.hat' matrix from equation 16
# well, rather semi-analytic: we compute the jacobian of the analytic
# Delta matrix
coef2JacobianSigma2 <- function(x) {
    lavmodel <- lav_model_set_parameters(fit@Model, x = x)
    DeltaS <- lavaan:::computeDelta(lavmodel)[[1]]
    lav_matrix_vec(t(DeltaS)) # note the transpose!
}
H.hat <- numDeriv:::jacobian(func = coef2JacobianSigma2, x = coef(fit))
# note: 1620 rows!!

# first derivative d F.ML / d c(Mu, vech(Sigma))
D.h1 <- -1/N * as.matrix(lavaan:::lav_mvnorm_dlogl_dmu_dvechSigma(Y = as.matrix(Data), Mu = h0.Mu, Sigma = h0.Sigma))

# note: t(D.h1) %*% Delta is the gradient wrt theta
D.theta <- as.matrix(lavaan:::lav_model_gradient(lavmodel = fit@Model, lavsamplestats = fit@SampleStats, lavdata = fit@Data))
tmp <- t(D.h1) %*% Delta
t(tmp) - D.theta

# identity matrix of size q = length(coef(fit))
Iq <- diag(length(coef(fit)))

# the second term in eq 15:
term2 <- (t(D.h1) %x% Iq) %*% H.hat

# analyic Hessian using eq 15 in Savalei & Rosseel (2022):
# note that we use '+' here (because we use F.ML, not the loglikelihood, there
# is a sign switch)
H.analytic <- A.theta.h1 + term2

# check
range(A.theta.num - H.analytic)
# -2.015421e-09  2.695497e-09

