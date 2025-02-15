test_that("Growth curve rel = 0.5", {
  the_seed <- 1234
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

  dist <- "Normal"
  nobs <- 15
  set.seed(the_seed)

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


  }

  data <- as.data.frame(t(lambda %*% factors + e))
  colnames(data) <- paste0("Day", 0:9)

  dat <- gen_data_growth(n = nobs, rel = 0.5, dist = dist, seed = the_seed)
  attr(dat, "truth") <- NULL
  attr(dat, "dist") <- NULL

  expect_equal(data, dat, ignore_attr = FALSE)
})

test_that("Growth curve rel = 0.8", {
  the_seed <- 1234
  dist <- "Kurtosis"
  nobs <- 20

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

  set.seed(the_seed)

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


  }

  data <- as.data.frame(t(lambda %*% factors + e))
  colnames(data) <- paste0("Day", 0:9)

  dat <- gen_data_growth(n = nobs, rel = 0.8, dist = dist, seed = the_seed)
  attr(dat, "truth") <- NULL
  attr(dat, "dist") <- NULL

  expect_equal(data, dat, ignore_attr = FALSE)
})

test_that("Two factor SEM rel = 0.5", {

  the_seed <- 789
  dist <- "NonNormal"
  nobs <- 50

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
  diag(reliability) <- rel <- 0.5

  # Define factor loadings, factor variances and residual variances
  lambda[1:6, 1] <- c(1, 0.7, 0.6, 0, 0, 0)
  lambda[1:6, 2] <- c(0, 0, 0, 1, 0.7, 0.6)
  diag(psi) <- c(1, 1)
  beta[2, 1] <- 0.25
  diag(theta) <- diag((lambda %*% psi %*% t(lambda)) %*%
                        solve(reliability) - (lambda %*% psi %*% t(lambda)))


  set.seed(the_seed)

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

  }

  data <- as.data.frame(t(lambda %*% factors + e))
  colnames(data) <- c("x1", "x2", "x3", "y1", "y2", "y3")

  dat <- gen_data_twofac(n = nobs, rel = rel, dist = "Non-normal",
                         seed = the_seed)
  attr(dat, "truth") <- NULL
  attr(dat, "dist") <- NULL

  expect_equal(data, dat, ignore_attr = FALSE)
})

test_that("Two factor SEM rel = 0.8", {

  the_seed <- 789
  dist <- "Normal"
  nobs <- 100

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
  diag(reliability) <- rel <- 0.8

  # Define factor loadings, factor variances and residual variances
  lambda[1:6, 1] <- c(1, 0.7, 0.6, 0, 0, 0)
  lambda[1:6, 2] <- c(0, 0, 0, 1, 0.7, 0.6)
  diag(psi) <- c(1, 1)
  beta[2, 1] <- 0.25
  diag(theta) <- diag((lambda %*% psi %*% t(lambda)) %*%
                        solve(reliability) - (lambda %*% psi %*% t(lambda)))

  set.seed(the_seed)

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

  }

  data <- as.data.frame(t(lambda %*% factors + e))
  colnames(data) <- c("x1", "x2", "x3", "y1", "y2", "y3")

  dat <- gen_data_twofac(n = nobs, rel = rel, dist = "Normal", seed = the_seed)
  attr(dat, "truth") <- NULL
  attr(dat, "dist") <- NULL

  expect_equal(data, dat, ignore_attr = FALSE)
})
