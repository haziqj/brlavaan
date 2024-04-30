# This script contains the code for the pairwise likelihood function.

library(tidyverse)
library(mvnfast)
library(mnormt)

# Function to generate data ----------------------------------------------------
gen_data_bin <- function(n = 1000, seed = NULL) {
  set.seed(seed)

  # Set up the loadings and covariance matrices --------------------------------
  Lambda      <- matrix(c(0.80, 0.70, 0.47, 0.38, 0.34), ncol = 1)
  neta        <- ncol(Lambda)  # q
  nitems      <- nrow(Lambda)  # p
  Psi         <- 1
  Theta       <- matrix(0, nrow = nitems, ncol = nitems)
  diag(Theta) <- 1 - diag(Lambda %*% Psi %*% t(Lambda)) # nolint
  tau         <- c(-1.43, -0.55, -0.13, -0.72, -1.13)

  # Generate the data ----------------------------------------------------------
  eta     <- mvnfast::rmvn(n = n, mu = rep(0, neta), sigma = Psi)
  epsilon <- mvnfast::rmvn(n = n, mu = rep(0, nitems), sigma = Theta)
  ystar <- tcrossprod(eta, Lambda) + epsilon

  y <-
    {1 * (ystar > matrix(tau, nrow = n, ncol = nitems, byrow = TRUE))} |>
    as.data.frame() |>
    mutate(across(everything(), \(x) ordered(x, levels = c(0, 1))))
  colnames(y) <- paste0("y", seq_len(nitems))

  as_tibble(y)
}

# Function to compute log-likelihood value -------------------------------------
create_pairwise_table <- function(data, wt = NULL) {

  p <- ncol(data)

  if (is.null(wt)) {
    data$wt <- 1
  } else {
    data$wt <- wt
  }

  n_combinations <- choose(p, 2)
  result_list <- vector("list", length = n_combinations)

  # Pre-compute column names
  col_names <- paste0("y", 1:p)

  k <- 1
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      var1 <- col_names[i]
      var2 <- col_names[j]

      # Create aggregate data frame
      agg_data <- data[c(var1, var2, "wt")]
      names(agg_data) <- c("yi", "yj", "freq")

      # Using tapply for aggregation
      freq <- tapply(agg_data$freq, list(agg_data$yi, agg_data$yj), length)
      wtfreq <- tapply(agg_data$freq, list(agg_data$yi, agg_data$yj), sum)

      # Convert to data frame and calculate proportions
      freq_df <- as.data.frame(as.table(freq))
      freq_df$wtfreq <- as.vector(wtfreq)
      freq_df$prop <- freq_df$Freq / sum(freq_df$Freq)
      freq_df$wtprop <- freq_df$wtfreq / sum(freq_df$wtfreq)

      # Store the result
      result_list[[k]] <- list(i = i, j = j, table = freq_df)
      k <- k + 1
    }
  }

  # Convert result list to a data frame or other desired format
  out <- do.call(rbind, lapply(result_list, as.data.frame))
  names(out) <- c("i", "j", "yi", "yj", "freq", "wtfreq", "prop", "wtprop")
  out$pattern <- paste0(out$yi, out$yj)
  out
}

# Pairwise likelihood function -------------------------------------------------
pl_fn <- function(theta, data, wt = NULL) {

  # 1. Create pairwise table x
  # 2. Get Var(ystar) x
  # 3. Mutate pattern probabilities
  # 4. Output the pairwise likelihood

  # Extract theta --------------------------------------------------------------
  lambdas <- matrix(theta[1:5], ncol = 1)  # the first 5 elements are loadings
  tau     <- theta[6:10] # the next 5 elements are thresholds
  Psi     <- 1  # no correlations in this 1-factor model

  # Compute Var(ystar) ---------------------------------------------------------
  neta <- ncol(Lambda) # nolint
  nitems <- nrow(Lambda)
  Theta <- matrix(0, nrow = nitems, ncol = nitems)
  diag(Theta) <- 1 - diag(Lambda %*% Psi %*% t(Lambda))
  Vy <- Lambda %*% Psi %*% t(Lambda) + Theta

  # Model probabilities --------------------------------------------------------
  calc_model_pairwise_prob <- function(i, j, pattern) {
    Vy_small <- Vy[c(i, j), c(i, j)]

    if (pattern == "00") {
      lower <- c(-Inf, -Inf)
      upper <- c(tau[i], tau[j])
    }
    if (pattern == "01") {
      lower <- c(-Inf, tau[j])
      upper <- c(tau[i], Inf)
    }
    if (pattern == "10") {
      lower <- c(tau[i], -Inf)
      upper <- c(Inf, tau[j])
    }
    if (pattern == "11") {
      lower <- c(tau[i], tau[j])
      upper <- c(Inf, Inf)
    }

    mnormt::sadmvn(lower = lower, upper = upper, mean = rep(0, 2),
                   varcov = Vy_small)
  }
  res <- create_pairwise_table(data, wt)
  res$prob <- with(res, mapply(calc_model_pairwise_prob, i, j, pattern))

  # Pairwise log-likelihood ----------------------------------------------------
  # sum(res$wtfreq * log(res$prob))
  with(res, sum(wtfreq * log(wtprop / prob)))

}