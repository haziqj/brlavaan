library(lavaan)
library(tidyverse)
library(furrr)
plan(multisession, workers = parallel::detectCores() - 2)
source("R/factor_models/10-ebrm_explicit.R")
source("R/factor_models/20-jackknife_bootstrap.R")
source("R/factor_models/33-pop_param.R")
B <- 1000  # no of sims
NB <- 500  # bootstrap samples

growth_sim_fn <- function(i = 1, n = 15, ver = 1, nboot = 5) {
  # Generate and fit data
  dat <- simulateData(model = txt_mod_growth(ver), sample.nobs = n)
  fit <- sem(
    model = "
    i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4 + 1*y5 + 1*y6 + 1*y7 + 1*y8 + 1*y9 + 1*y10
    s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4 + 4*y5 + 5*y6 + 6*y7 + 7*y8 + 8*y9 + 9*y10
    ",
    data = dat
  )

  # Bias reducing methods
  .jack <- rb_jack(fit, dat)
  if(inherits(.jack, "try-error")) jack_time <- difftime(NA, Sys.time())
  else jack_time <- attr(.jack, "timing")

  .boot <- try(rb_boot(fit, dat, nboot), silent = TRUE)
  if(inherits(.boot, "try-error")) boot_time <- difftime(NA, Sys.time())
  else boot_time <- attr(.boot, "timing")

  .ebrm <- try(rb_empr(fit))
  if(inherits(.ebrm, "try-error")) ebrm_time <- difftime(NA, Sys.time())
  else ebrm_time <- attr(.ebrm, "timing")

  tibble(
    i     = i,
    n     = n,
    rel   = ifelse(ver == 1, "0.8", "0.5"),
    par   = names(coef(fit)),
    truth = truth_growth(ver),
    ml    = as.numeric(coef(fit)),
    jack  = as.numeric(.jack),
    boot  = as.numeric(.boot),
    ebrm  = as.numeric(.ebrm),
    # Convergences
    conv_ml   = fit@Fit@converged,
    conv_jack = sum(attr(.jack, "meta")$ok),
    conv_boot = sum(attr(.boot, "meta")$ok),
    conv_ebrm = !inherits(.ebrm, "try-error"),
    # Timings
    time_jack = jack_time,
    time_boot = boot_time,
    time_ebrm = ebrm_time
  )
}

res_growth <- list()
i <- 1
for (samp_size in c(15, 20, 50, 100, 1000)) {
  for (ver in 1:2) {
    rel <- ifelse(ver == 1, "0.8", "0.5")
    cli::cli_alert_info("Running growth curve sims (rel = {rel} / n = {samp_size})")

    res_growth[[i]] <-
      furrr::future_map(
        .x = 1:B,
        .f = \(x) growth_sim_fn(i = x, n = samp_size, ver = ver, nboot = NB),
        .progress = TRUE,
        .options = furrr_options(seed = TRUE)
      )

    i <- i + 1
    cat("\n")
  }
}

save(res_growth, file = "R/factor_models/sim_growth_curve.RData")

res_growth_all <-
  lapply(res_growth, \(x) do.call(rbind, x)) |>
  do.call(what = rbind) |>
  mutate(par2 = case_when(
    par == "i~~i" ~ "Psi11",
    par == "s~~s" ~ "Psi22",
    par == "i~~s" ~ "Psi12",
    TRUE ~ "theta"
  )) |>
  mutate(
    par2 = factor(par2, levels = c("theta", "Psi11", "Psi22", "Psi12"))
  )

# Convergence failures and timing
res_growth_all |>
  summarise(
    across(starts_with("conv"), first),
    across(starts_with("time"), first),
    .by = c(i, rel, n)
  ) |>
  summarise(
    fail_ml = 1 - mean(conv_ml),
    fail_jack = 1 - mean(conv_jack / n),
    fail_boot = 1 - mean(conv_boot / NB),
    across(starts_with("time"), \(x) mean(x, na.rm = TRUE)),
    .by = c(rel, n)
  ) |>
  arrange(desc(rel), n)

# Relative mean bias table
res_growth_summary <-
  res_growth_all |>
  summarise(
    across(c(ml, ebrm, jack, boot), \(x) mean((x - truth)/truth, na.rm = TRUE)),
    .by = c(rel, n, par2)
  ) |>
  pivot_longer(
    cols = ml:boot,
    names_to = "method",
    values_to = "rel_bias"
  ) |>
  mutate(method = factor(method, c("ml", "jack", "boot", "ebrm"))) |>
  arrange(desc(rel), n, par2, method)

res_growth_summary |>
  pivot_wider(
    names_from = n,
    names_prefix = "n = ",
    values_from = rel_bias
  )

# Relative mean bias plot
res_growth_summary |>
  mutate(
    x = as.numeric(as.factor(n)),
    rel = paste0("Reliability = ", rel)
  ) |>
  ggplot(aes(x, rel_bias, color = method, shape = method)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "gray") +
  geom_point() +
  geom_line() +
  scale_x_continuous(
    breaks = c(1, 2, 3, 4, 5),
    labels = c(15, 20, 50, 100, 1000),
    name = "Sample size"
  ) +
  scale_y_continuous(labels = scales::percent, name = "Relative mean bias") +
  ggsci::scale_colour_npg(name = NULL) +
  scale_shape_manual(values = c(15, 5, 9, 16), name = NULL) +
  facet_grid(par2 ~ rel) +
  theme_bw() +
  theme(
    legend.position = "top"
  )
