library(tidyverse)
library(lavaan.bingof)
library(kableExtra)

## A dummy fit to get par names
par_names <- names(theta_hat)

## Simulation studies
set.seed(123)
n_simu <- 100
samp_size <- 1000
n_cores <- parallel::detectCores() - 2
model_no <- 1
theta_true <- c(lavaan.bingof:::get_Lambda(model_no),
                lavaan.bingof:::get_tau(model_no))
mod <- txt_mod(model_no)
datasets <- replicate(n_simu, gen_data_bin(model_no, n = samp_size),
                      simplify = FALSE)

failed <- function(j, n_simu) cat("Sample:", j, "/", n_simu, ": Failed\n")
done <- function(j, n_simu) cat("Sample:", j, "/", n_simu, ": Done\n")

results <- parallel::mclapply(seq.int(n_simu), function(j) {
  fit_lav <- try(sem(mod, datasets[[j]], std.lv = TRUE, estimator = "PML"),
                 silent = TRUE)
  if (is(fit_lav, "try-error")) {
    mpl <- eRBM <- rep(NA, length(theta_true))
    failed(j, n_simu)
  } else {
    mpl <- coef(fit_lav)
    A <- try(AAA(mpl), silent = TRUE)
    Hinv <- try(solve(HHH(mpl, unit = FALSE)), silent = TRUE)
    if (is(A, "try-error") | is(Hinv, "try-error")) {
      eRBM <- rep(NA, length(theta_true))
      failed(j, n_simu)
    } else {
      eRBM <- drop(mpl + Hinv %*% A)
      done(j, n_simu)
    }
  }
  data.frame(estimate = c(mpl, eRBM),
             method = rep(c("MCL", "eRBM"), each = length(par_names)),
             truth = c(theta_true, theta_true),
             sample = j,
             parameter = rep(par_names, 2))
}, mc.cores = n_cores)



res <- do.call("rbind", results)
res |>
  group_by(parameter, method) |>
  filter(!is.na(estimate)) |>
  summarize(bias = mean(estimate - truth),
            sd = sd(estimate),
            abs_bias = mean(abs(estimate - truth)),
            rmse = sqrt(mean((estimate - truth)^2)),
            min = min(estimate),
            max = max(estimate),
            n_used = n())
  # kbl(booktabs = TRUE, digits = 2, linesep = "")
