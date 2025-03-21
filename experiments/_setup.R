library(brlavaan)
library(tidyverse)
library(furrr)
theme_set(theme_bw())
load(here::here("experiments/simu_id.RData"))

ncores <- future::availableCores() - 1
future::plan(multisession, workers = ncores)
B <- 500  # Number of simulations

## ----- Simulation function ---------------------------------------------------
sim_fun <- function(
    dist = c("Normal", "Kurtosis", "Non-normal"),
    model = c("growth", "twofac"),
    rel = 0.8,
    n = 15,
    nsimu = 1,
    nboot = 500L,
    lavsim = FALSE,
    whichsims = c("ML", "eRBM", "iRBM", "Ozenne", "REML"),
    bounds = "standard",
    data_scale = 1,
    seeds = NULL
  ) {
  dist <- match.arg(dist)
  model <- match.arg(model)
  rel <- match.arg(as.character(rel), c("0.8", "0.5"))
  whichsims <- match.arg(whichsims, several.ok = TRUE,
                         c("ML", "eRBM", "iRBM", "Ozenne", "REML", "JB", "BB"))

  # Check seeds
  if (!is.null(seeds)) {
    if (length(seeds) != nsimu) {
      cli::cli_abort("Length of seeds must be equal to nsimu")
    }
  }

  # Generate all data ----------------------------------------------------------
  if (model == "growth") {
    gen_data <- gen_data_growth
    txt_mod <- txt_mod_growth
    lavfun <- "growth"
  } else if (model == "twofac") {
    gen_data <- gen_data_twofac
    txt_mod <- txt_mod_twofac  # this has arg meanstructure = FALSE by default
    lavfun <- "sem"
  } else {
    cli::cli_abort("Unknown model: {model}")
  }
  datasets <- NULL
  if (!is.null(seeds)) {
    datasets <- purrr::map(
      .x = seeds,
      .f = \(s) gen_data(n = n, rel = rel, dist = dist, lavsim = lavsim,
                         scale = data_scale, seed = s)
    )
  }
  mod <- txt_mod(rel)

  # Single run function --------------------------------------------------------
  single_sim <- function(j) {
    if (is.null(datasets)) {
      dat <- gen_data(n = n, rel = rel, dist = dist, lavsim = lavsim,
                      scale = data_scale)
    } else {
      dat <- datasets[[j]]
    }
    true_vals <- truth(dat)

    fitsemargs <- list(
      model = mod,
      data = dat,
      rbm = "none",
      plugin_pen = NULL,
      debug = FALSE,
      lavfun = lavfun,
      maxgrad = FALSE,
      fn.scale = 1,
      bounds = bounds,
      start = true_vals
    )

    nsimtypes <- length(whichsims) + 1  # to account for lav below, which will always be fit
    fit_list <- list()

    # OUR SIMS -----------------------------------------------------------------
    if ("ML" %in% whichsims) {
      fit_list$ML <- try(do.call(fit_sem, fitsemargs))
    }
    if ("eRBM" %in% whichsims) {
      fitsemargs$rbm <- "explicit"
      fit_list$eRBM <- try(do.call(fit_sem, fitsemargs))
    }
    if ("iRBM" %in% whichsims) {
      fitsemargs$rbm <- "implicit"
      fit_list$iRBM <- try(do.call(fit_sem, fitsemargs))
    }

    # D&R SIMS -----------------------------------------------------------------
    if (model == "growth") {
      fit_lav <- try(growth(mod, dat, start = true_vals, bounds = bounds, ceq.simple = TRUE))
    } else if (model == "twofac") {
      fit_lav <- try(sem(mod, dat, start = true_vals, bounds = bounds))
    }
    if (!inherits(fit_lav, "try-error")) {
      out <- list()
      out$coefficients <- as.numeric(coef(fit_lav))
      out$stderr <- partable(fit_lav)$se[partable(fit_lav)$free > 0]
      out$timing <- fit_lav@timing$total
      out$converged <- fit_lav@optim$converged
      out$Sigma <- vcov(fit_lav)
      fit_list$lav <- out
    } else {
      fit_list$lav <- NA
    }
    lavok <- ifelse(length(fit_list$lav) == 1, !is.na(fit_list$lav), TRUE)

    if ("Ozenne" %in% whichsims) {
      if (lavok) {
        start_time <- proc.time()
        out <- try(brlavaan:::yr_ozenne_bcs_growth(fit_lav, return.se = TRUE),
                   silent = TRUE)
        elapsed_time <- proc.time() - start_time
        elapsed_time <- elapsed_time["elapsed"]
      } else{
        out <- try(log("a"), silent = TRUE)
      }
      if (!inherits(out, "try-error")) {
        names(out) <- c("coefficients", "stderr")
        out$timing <- elapsed_time
        out$converged <- fit_list$lav$converged
        out$Sigma <- brlavaan:::get_Sigma(out$coefficients, fit_lav)
        fit_list$Ozenne <- out
      } else {
        fit_list$Ozenne <- NA
      }
    }
    if ("JB" %in% whichsims) {
      if (lavok) {
        start_time <- proc.time()
        jack_out <- try(brlavaan:::sd_jack_bc(fit_lav), silent = TRUE)
        elapsed_time <- proc.time() - start_time
        elapsed_time <- elapsed_time["elapsed"]
      } else{
        jack_out <- try(log("a"), silent = TRUE)
      }
      if (!inherits(jack_out, "try-error")) {
        out <- list()
        out$coefficients <- as.numeric(jack_out$theta.bc)
        out$stderr <- jack_out$jack.se
        out$timing <- elapsed_time
        out$converged <- fit_list$lav$converged
        out$Sigma <- brlavaan:::get_Sigma(out$coefficients, fit_lav)
        fit_list$JB <- out
      } else {
        fit_list$JB <- NA
      }
    }
    if ("BB" %in% whichsims) {
      if (lavok) {
        start_time <- proc.time()
        boot_out <- try(brlavaan:::sd_boot_bc(fit_lav, ndraws = nboot,
                                   jack.resamples = jack_out$jack.resamples),
                        silent = TRUE)
        elapsed_time <- proc.time() - start_time
        elapsed_time <- elapsed_time["elapsed"]
      } else{
        boot_out <- try(log("a"), silent = TRUE)
      }
      if (!inherits(boot_out, "try-error")) {
        out <- list()
        out$coefficients <- as.numeric(boot_out$theta.bc)
        out$stderr <- boot_out$boot.se
        out$timing <- elapsed_time
        out$converged <- fit_list$lav$converged
        out$Sigma <- brlavaan:::get_Sigma(out$coefficients, fit_lav)
        fit_list$BB <- out
      } else {
        fit_list$BB <- NA
      }
    }
    if ("REML" %in% whichsims & model == "growth") {
      start_time <- proc.time()
      out <- try(brlavaan:::dr_reml_growth(dat), silent = TRUE)
      elapsed_time <- proc.time() - start_time
      elapsed_time <- elapsed_time["elapsed"]
      if (!inherits(out, "try-error")) {
        names(out) <- c("coefficients", "stderr")
        out$timing <- elapsed_time
        out$converged <- fit_list$ML$converged
        out$Sigma <- brlavaan:::get_Sigma(out$coefficients, fit_lav)
        fit_list$REML <- out
      } else {
        fit_list$REML <- NA
      }
    }

    tibble::tibble(
      seed = seeds[j],
      sim = j,
      dist = dist,
      model = model,
      rel = rel,
      n = n,
      method = names(fit_list),
      est = lapply(fit_list, purrr::possibly(\(x) x$coefficients, NA_real_)),
      se = lapply(fit_list, purrr::possibly(\(x) x$stderr, NA_real_)),
      truth = rep(list(true_vals), nsimtypes),
      timing = sapply(fit_list, purrr::possibly(\(x) x$timing, NA_real_)),
      converged = sapply(fit_list, purrr::possibly(\(x) x$converged, NA)),
      Sigma_OK = sapply(fit_list, purrr::possibly(\(x) !brlavaan:::check_mat(x$Sigma), NA))
    )
  }

  # Run simulation -------------------------------------------------------------
  simu_res <- furrr::future_map(
    seq_len(nsimu),
    purrr::possibly(single_sim),
    .progress = TRUE,
    .options = furrr::furrr_options(seed = TRUE)
  )
  # simu_res
  do.call("rbind", simu_res)
}
