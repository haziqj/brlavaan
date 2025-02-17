##==============================================================================
## Project:     Resampling based bias correction for small sample SEM
##
## Script:      (Simplified) implementation of the bias correction method
##              described in Ozenne et al., 2020; only for models without
##              subject-dependent X_i covariates.
##==============================================================================

yr_ozenne_bcs_twofac <- function(lavobject,
                          verbose = FALSE,
                          max.iter = 100L,
                          tol = 1e-05,
                          return.se = FALSE) {

  # meanstructure = TRUE is needed (for Delta.mu)
  stopifnot(lavobject@Model@meanstructure == TRUE,
            lavobject@Model@categorical == FALSE,
            lavobject@Model@conditional.x == FALSE)
  stopifnot(lavobject@Data@ngroups == 1L,
            lavobject@Data@nlevels == 1L)

  lavsamplestats = lavobject@SampleStats
  lavoptions     = lavobject@Options
  lavdata        = lavobject@Data
  lavpartable    = lavobject@ParTable
  lavmodel       = lavobject@Model

  conditional.x  = lavmodel@conditional.x

  # create selection matrix (assuming a single group)
  xm2select <- function(lavmodel, target = "theta") {
      target.idx <- which(names(lavmodel@GLIST) == target)[1] # first group
      # dimension
      P <- nrow(lavmodel@GLIST[[target.idx]])
      unique.idx <- unique(lavmodel@x.free.idx[[target.idx]])
      row.idx <- match(lavmodel@x.free.idx[[target.idx]], unique.idx)
      out <- matrix(0, nrow = P*P, ncol = length(unique.idx))
      IDX <- cbind(lavmodel@m.free.idx[[target.idx]], row.idx)
      out[IDX] <- 1
      out
  }

  # L.PSI, L.THETA and G2
  MLIST <- lavmodel@GLIST
  nvar <- nrow(MLIST$lambda)
  nfac <- ncol(MLIST$lambda)
  N <- lavaan::lavTech(lavobject, "nobs")
  L.PSI <- xm2select(lavmodel, target = "psi")
  L.THETA <- xm2select(lavmodel, target = "theta")
  G2 <- lavaan::lav_matrix_duplication_ginv_pre(L.THETA)
  if(!is.null(MLIST$beta)) {
    IB.inv <- lavaan___.internal_get_IB.inv(MLIST = MLIST)
    Lambda.IB.inv <- MLIST$lambda %*% IB.inv
    LxL <- Lambda.IB.inv %x% Lambda.IB.inv
  } else {
    LxL <- MLIST$lambda %x% MLIST$lambda
  }
  G1 <- lavaan::lav_matrix_duplication_ginv_pre(LxL) %*% L.PSI
  G <- cbind(G1, G2)

  if(conditional.x) {
    Sigma.k <- Sigma <- lavaan::lavTech(lavobject, "implied")[[1]]$res.cov
  } else {
    Sigma.k <- Sigma <- lavaan::lavTech(lavobject, "implied")[[1]]$cov
  }

  lavmodel.k <- lavmodel
  VCOV.k <- VCOV <- lavaan::lavTech(lavobject, "vcov")
  Delta.mu.k <- Delta.mu <- lav_model_dmu_dtheta(lavobject@Model)
  old.k <- lavaan::lav_matrix_vech(Sigma.k)

  for(k in seq_len(max.iter)) {

      if(verbose) {
          cat("iteration = ", k)
      }

      # a. + b. compute the bias for each subject, take the average
      # here, we assume every 'i' has the same value for Delta.mu (for now)
      Sigma.bias <- Delta.mu.k %*% VCOV.k %*% t(Delta.mu.k) / N

      # c. update Sigma
      Sigma.k <- Sigma + Sigma.bias

      # check for convergence
      new.k <- lavaan::lav_matrix_vech(Sigma.k)
      rms.k <- sqrt(mean((old.k - new.k) * (old.k - new.k)))
      if(rms.k < tol) {
          if(verbose) {
              cat(" ... rms.k = ", rms.k, " -- converged.\n")
          }
      } else {
          old.k <- new.k
          if(verbose) {
              cat(" ... rms.k = ", rms.k, "\n")
          }
      }

      # d. find theta.var.k using OLS
      out <- solve(crossprod(G)) %*% t(G) %*% lavaan::lav_matrix_vech(Sigma.k)

      # reconstruct psi
      MLIST$psi <- matrix(L.PSI %*% out[1:ncol(L.PSI)], nfac, nfac)

      # reconstruct theta
      MLIST$theta <- matrix(L.THETA %*% out[(1+ncol(L.PSI)):length(out)],
                            nvar, nvar)

      # updated parameter vector and lavmodel
      theta.updated <- lavaan::lav_model_get_parameters(lavmodel,
                                                GLIST = MLIST,
                                                type = "free")

      # recreate lavmodel with update parameter vector
      lavmodel.k <- lavaan::lav_model_set_parameters(lavmodel, x = theta.updated)

      # e. + f. update VCOV.k
      VCOV.k <- lavaan___lav_model_vcov(lavmodel = lavmodel.k,
                                        lavsamplestats = lavsamplestats,
                                        lavoptions = lavoptions,
                                        lavdata = lavdata) * N

      # update Delta.mu (needed?)
      Delta.mu.k <- lavaan___computeDelta(lavmodel.k)[[1]][1:nvar, ,
                                                           drop = FALSE]

    } # iterations

  cor.se <- ifelse(test = grepl("~1", names(lavaan::coef(lavobject))),
                   yes = sqrt(diag(VCOV.k)/N),
                   no = sqrt(diag(VCOV.k)/(N-1)))

  if(return.se) {
    list("est.corrected" = theta.updated,
         "se" = cor.se)
  } else {
    return("est.corrected" = theta.updated)

  }
}


# compute 'Delta.mu', even if meanstructure = FALSE
# means in the rows, theta in the columns
lav_model_dmu_dtheta <- function(lavmodel) {

  GLIST <- lavmodel@GLIST

  NCOL <- lavmodel@nx.free
  m.el.idx <- x.el.idx <- vector("list", length=length(GLIST))
  for(mm in 1:length(GLIST)) {
      m.el.idx[[mm]] <- lavmodel@m.free.idx[[mm]]
      x.el.idx[[mm]] <- lavmodel@x.free.idx[[mm]]
      # handle symmetric matrices
      if(lavmodel@isSymmetric[mm]) {
          # since we use 'x.free.idx', only symmetric elements
          # are duplicated (not the equal ones, only in x.free.free)
          dix <- duplicated(x.el.idx[[mm]])
          if(any(dix)) {
              m.el.idx[[mm]] <- m.el.idx[[mm]][!dix]
              x.el.idx[[mm]] <- x.el.idx[[mm]][!dix]
          }
      }
  }

  Delta.mu <- matrix(0, nrow = lavmodel@nvar[1], ncol = NCOL)

  # single group (for now)!
  for(mm in seq_len(lavmodel@nmat)) {
      mname <- names(lavmodel@GLIST)[mm]

      # skip empty ones
      if(!length(m.el.idx[[mm]])) next

      Delta.mu[, x.el.idx[[mm]]] <-
          lavaan___derivative.mu.LISREL(m = mname,
                                        idx = m.el.idx[[mm]],
                                        MLIST = GLIST)
    } # mm

    Delta.mu
}
