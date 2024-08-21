library(tidyverse)
theme_set(theme_bw())
library(lavaan)
library(furrr)
ncores <- parallel::detectCores() - 2
future::plan(multisession, workers = ncores)
source("R/41-implicit_wrap.R")

theta_true <- c(
  0.8, 0.6, 0.8, 0.6,  # lambda
  0.3,  # beta
  rep(1, 6),  # theta_eps
  rep(1, 2)  # psi
)
m <- length(theta_true)

B <- 10  # no of sims
res <- list()
nvec <- c(15, 20, 50, 100, 1000)

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
  res[[i]] <- future_map(
    datlist,
    imp_twofac,
    .progress = TRUE
  )
  cat("\n")
}
names(res) <- nvec

res_summ <-
  lapply(
    res,
    \(x) {
      ress <-
        as_tibble(do.call("rbind", x)) |>
        mutate(
          iRBMp_orig = iRBMp,
          ok = sapply(iRBMp, \(x) !inherits(x, "try-error")),
          iRBMp = sapply(iRBMp, \(x) ifelse(inherits(x, "try-error"), NA, list(x))),
          nok = sum(ok),
          truth = list(theta_true)
        ) |>
        filter(ok)

      if (nrow(ress) > 0) {
        ress |>
          unnest(c(ml, iRBMp, truth)) |>
          mutate(param = names(ml)) |>
          select(param, truth, ml, iRBMp, nok) |>
          mutate(
            across(c(ml, iRBMp), \(x) x - truth),
            type = case_when(
              grepl("=~", param) ~ "loadings",
              grepl("~~", param) ~ "variance",
              TRUE ~ "regression"
            )
          )
      } else {
        ress
      }
    }
  )
res_summ <- imap(res_summ, ~ mutate(.x, n = .y))
res_summ <- do.call("rbind", res_summ)
res_summ$n <- factor(res_summ$n, levels = nvec)

# Results ----------------------------------------------------------------------
res_summ |>
  pivot_longer(c(ml, iRBMp), names_to = "method", values_to = "bias") |>
  ggplot(aes(as.numeric(n), abs(bias), col = method)) +
  geom_point(size = 0.5, alpha = 0.5, position = position_dodge(width = 0.1)) +
  geom_smooth(se = FALSE) +
  scale_y_log10() +
  scale_x_continuous(breaks = 1:5, labels = nvec) +
  facet_grid(type ~ .) +
  labs(
    x = "Sample size",
    y = "Absolute bias",
    col = "Method"
  )

# Mean absolute bias
res_summ |>
  group_by(n) |>
  summarise(
    across(c(ml, iRBMp), ~ mean(abs(.x), na.rm = TRUE)),
    .groups = "drop"
  )



