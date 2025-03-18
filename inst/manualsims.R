library(tidyverse)
library(furrr)
plan(multisession)

res <- future_map(1:100, \(i) {
  require(brlavaan)
  dat <- gen_data_growth(n = 25, rel = 0.8, dist = "Normal", scale = 1 / 10)
  mod <- txt_mod_growth(0.8)
  tru <- truth(dat)

  fit <- list()
  fit$ML <- growth(mod, dat, start = tru)
  fit$eBR <- fit_sem(mod, dat, "explicit", lavfun = "growth", start = tru)
  fit$iBR <- fit_sem(mod, dat, "implicit", lavfun = "growth", start = tru)
  fit$MLb  <- growth(mod, dat, start = tru, bounds = "standard")
  fit$eBRb <- fit_sem(mod, dat, "explicit", lavfun = "growth", start = tru, bounds = "standard")
  fit$iBRb <- fit_sem(mod, dat, "implicit", lavfun = "growth", start = tru, bounds = "standard")
  fit$eBRman <- fit_growth(mod, dat, "explicit", start = tru[1:6])
  fit$iBRman <- fit_growth(mod, dat, "implicit", start = tru[1:6])

  c(list(truth = tru[1:6]), lapply(fit, \(x) {
    conv <- if (inherits(x, "lavaan")) x@optim$converged else x$converged
    if (conv) coef(x)[1:6] else NA
  })) |>
    as.data.frame()
}, .progress = TRUE, .options = furrr_options(seed = TRUE))

map(res, \(x) rownames_to_column(x, "param")) |>
  bind_rows(.id = "sim") |>
  pivot_longer(cols = ML:iBRman, names_to = "method", values_to = "est") |>
  mutate(bias = est - truth) |>
  # filter(abs(bias) < 100) |>
  summarise(
    bias = mean(bias, na.rm = TRUE, trim = 0.05),
    .by = c(param, method)
  ) |>
  mutate(
    method = factor(method, levels = c("MLb", "eBRb", "iBRb")),
    param = factor(param, levels = c("i~~i", "i~1", "s~~s", "s~1", "i~~s", "v"))
  ) |>
  ggplot(aes(param, abs(bias), fill = method)) +
  geom_col(position = "dodge", width = 0.5) +
  # scale_y_log10() +
  # coord_cartesian(ylim = c(0, 1)) +
  theme_minimal()

map(res, \(x) rownames_to_column(x, "param")) |>
  bind_rows(.id = "sim") |>
  mutate(across(ML:iBRman, \(x) !is.na(x))) |>
  summarise(across(ML:iBRman, sum), .by = param)
