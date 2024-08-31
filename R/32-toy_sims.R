library(tidyverse)
library(lavaan)
library(furrr)
theme_set(theme_bw())
source("R/10-sem_rbm_functions.R")
source("R/20-gen_data.R")
source("R/21-sim_functions.R")

ncores <- parallel::detectCores() - 1
future::plan(multisession, workers = ncores)
B <- 1000  # no of sims

# Parameters -------------------------------------------------------------------
true_vals <- c(
  0.8, 0.6, 0.8, 0.6,  # lambda
  0.3,  # beta
  rep(1, 6),  # theta_eps
  rep(1, 2)  # psi
)
names(true_vals) <- c(
  "eta1=~y2", "eta1=~y3", "eta2=~y5", "eta2=~y6",
  "eta2~eta1",
  "y1~~y1", "y2~~y2", "y3~~y3", "y4~~y4", "y5~~y5", "y6~~y6",
  "eta1~~eta1", "eta2~~eta2"
)

res <- list()
nvec <- c(15, 20, 50, 100, 1000)

# Sim function -----------------------------------------------------------------
toysim <- function(dat) {
  mod <- "
    eta1 =~ y1 + y2 + y3
    eta2 =~ y4 + y5 + y6
    eta2 ~ eta1
  "

  fit_ML    <- fit_sem(mod, dat, method = "ML")
  fit_eRBM  <- fit_sem(mod, dat, method = "eRBM")
  fit_iRBM  <- fit_sem(mod, dat, method = "iRBM")
  fit_iRBMp <- fit_sem(mod, dat, method = "iRBMp")

  data.frame(
    n = n,
    param = rep(names(true_vals), 4),
    est = c(fit_ML$coef, fit_eRBM$coef, fit_iRBM$coef, fit_iRBMp$coef),
    se = c(fit_ML$stderr, fit_eRBM$stderr, fit_iRBM$stderr, fit_iRBMp$stderr),
    truth = rep(true_vals, 4),
    method = rep(c("ML", "eRBM", "iRBM", "iRBMp"), each = length(true_vals)),
    time = rep(c(fit_ML$time, fit_eRBM$time, fit_iRBM$time, fit_iRBMp$time), each = length(true_vals))
  )
}

for (i in seq_along(nvec)) {
  n <- nvec[i]

  datlist <- replicate(B, {
    simulateData(
      model = "
      eta1 =~ 1*y1 + 0.8*y2 + 0.6*y3
      eta2 =~ 1*y4 + 0.8*y5 + 0.6*y6
      eta2 ~ 0.3*eta1
    ",
      sample.nobs = n
    )
  }, simplify = FALSE)

  cli::cli_inform("[{i} / {length(nvec)}] Now running n = {n}\n")
  res[[i]] <- future_map(datlist, safely(toysim, NA), .progress = TRUE)
  cat("\n")
}

# Clean up
simu_res <- lapply(res, function(simu_res) {
  where_error <- which(sapply(simu_res, \(x) !is.null(x$error)))
  error_list <- lapply(simu_res[where_error], \(x) x$error)
  names(error_list) <- where_error

  list(
    simu_res = do.call("rbind", lapply(simu_res, \(x) x$result)),
    errors = error_list
  )
})

results <-
  do.call(rbind, lapply(simu_res, \(x) x$simu_res)) |>
  as_tibble() |>
  mutate(
    n = factor(n),
    method = factor(method, levels = c("ML", "eRBM", "iRBM", "iRBMp")),
    bias = est - truth,
    cov  = abs(bias) < 1.96 * se,
    type = case_when(
      grepl("=~", param) ~ "Loadings",
      grepl("~~", param) ~ "Variances",
      TRUE ~ "Regressions"
    )
  )

# Analysis
tab <-
  results |>
  filter(abs(est) < 1000) |>
  group_by(n, method, type) |>
  summarise(
    bias = mean(bias),
    mse = mean(bias ^ 2),
    .groups = "drop"
  ) |>
  pivot_wider(id_cols = c(n, type), names_from = method, values_from = c(bias, mse))

cutoff <- 1000
p <-
  results |>
  filter(abs(est) < cutoff) |>
  ggplot(aes(as.numeric(n), abs(bias), col = method)) +
  geom_point(size = 0.5, alpha = 0.1, position = position_dodge(width = 0.1)) +
  geom_line(
    data = results |>
      filter(abs(est) < cutoff) |>
      group_by(n, method) |>
      summarise(
        bias = mean(bias),
        mse = mean(bias ^ 2),
        .groups = "drop"
      ),
    aes(y = abs(bias)),
    linewidth = 1,
    alpha = 0.8
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", \(x) 10 ^ x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
  ) +
  scale_x_continuous(labels = nvec) +
  labs(
    y = "(Absolute) Mean bias",
    x = "Sample size",
    col = NULL
  ) +
  theme(legend.position = "top")



save(tab, p, file = "toysims.RData")








