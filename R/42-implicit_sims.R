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

B <- 100  # no of sims
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
              grepl("=~", param) ~ "Loadings",
              grepl("~~", param) ~ "Variances",
              TRUE ~ "Regressions"
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

# Mean absolute bias table
tab_mab <-
  res_summ |>
  group_by(n) |>
  summarise(
    across(c(ml, iRBMp), ~ mean(abs(.x), na.rm = TRUE)),
    nok = mean(nok),
    .groups = "drop"
  )

tab_mab_by_type <-
  res_summ |>
  group_by(n, type) |>
  summarise(
    across(c(ml, iRBMp), ~ mean(abs(.x), na.rm = TRUE)),
    nok = mean(nok),
    .groups = "drop"
  )

# Bias plot
p_bias <-
  res_summ |>
  pivot_longer(c(ml, iRBMp), names_to = "method", values_to = "bias") |>
  ggplot(aes(as.numeric(n), abs(bias), col = method)) +
  geom_point(size = 0.5, alpha = 0.3, position = position_dodge(width = 0.1)) +
  geom_line(
    data = pivot_longer(
      tab_mab_by_type,
      c(ml, iRBMp),
      names_to = "method",
      values_to = "bias"
    ),
    linewidth = 1
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", \(x) 10 ^ x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
  ) +
  scale_x_continuous(breaks = 1:5, labels = nvec) +
  facet_grid(type ~ .) +
  labs(
    x = "Sample size",
    y = "Absolute bias",
    col = "Method"
  )

save(p_bias, tab_mab, tab_mab_by_type, file = "R/prelim_implicit_sims.RData")
