skip()
library(tidyverse)
here::i_am("tests/testthat/test-DRfit.R")
dr_file1 <- here::here("experiments/GCM_est_combined_final.RData")
dr_file2 <- here::here("experiments/2FSEM_est_combined_final.RData")
if (!file.exists(dr_file1))
  download.file("https://osf.io/vjq5m/download", destfile = dr_file1)
if (!file.exists(dr_file2))
  download.file("https://osf.io/cw5b7/download", destfile = dr_file2)

test_that("Growth curve model results are similar with D&R", {

  load(here::here("experiments/GCM_est_combined_final.RData"))

  # Sample a couple of results
  set.seed(123)
  DR_Res <-
    Results |>
    filter(method == "MLB") |>
    filter(seed == 1235, dist == "Normal", rel == "REL80", nobs == 15) |>
    # slice_sample(n = 10) |>
    select(rel, dist, method, nobs, seed, starts_with("est")) |>
    nest(est = starts_with("est")) |>
    mutate(est = map(est, \(x) {
      colnames(x) <- gsub("est_", "", colnames(x))
      x
    }))

  for (i in seq_len(nrow(DR_Res))) {
    n <- as.numeric(as.character(DR_Res$nobs[i]))
    rel <- as.numeric(sub("REL", "", DR_Res$rel[i])) / 100
    dist <- as.character(DR_Res$dist[i])
    dist <- "Normal"
    if (dist == "NonNormal") dist <- "Non-normal"
    seed <- DR_Res$seed[i]
    estDR <- unlist(DR_Res$est[i])

    dat <- gen_data_growth(n = n, rel = rel, dist = dist, seed = seed)
    mod <- txt_mod_growth(rel)
    fit <- growth(mod, dat, start = truth_growth(rel), bounds = "standard", meanstructure = TRUE)
    est <- coef(fit)
    class(est) <- "numeric"

    expect_equal(est, estDR, tolerance = 1e-2, ignore_attr = TRUE)
  }
})


test_that("Two factor results are similar with D&R", {

  load(here::here("experiments/2FSEM_est_combined_final.RData"))

  # Sample a couple of results
  DR_Res <-
    Results |>
    filter(method == "MLB") |>
    slice_sample(n = 10) |>
    select(rel, dist, method, nobs, seed, starts_with("est")) |>
    nest(est = starts_with("est")) |>
    mutate(est = map(est, \(x) {
      colnames(x) <- gsub("est_", "", colnames(x))
      x
    }))

  for (i in seq_len(nrow(DR_Res))) {
    n <- as.numeric(as.character(DR_Res$nobs[i]))
    rel <- as.numeric(sub("REL", "", DR_Res$rel[i])) / 100
    dist <- as.character(DR_Res$dist[i])
    if (dist == "NonNormal") dist <- "Non-normal"
    seed <- DR_Res$seed[i]

    set.seed(seed)
    dat <- gen_data_twofac(n = n, rel = rel, dist = dist)
    mod <- txt_mod_twofac(rel)
    fit <- sem(mod, dat, start = truth_twofac(rel, meanstructure = TRUE),
               bounds = "standard", meanstructure = TRUE)
    est <- coef(fit)
    class(est) <- "numeric"

    expect_equal(est, unlist(DR_Res$est[i]), tolerance = 1e-4)
  }
})
