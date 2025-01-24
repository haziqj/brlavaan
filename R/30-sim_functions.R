#' Conduct bias reduction simulation study
#'
#' Conduct bias reduction simulation study
#'
#' Suggest to enable parallel processing by calling
#' `future::plan("multisession", workers = future::availableCores() - 1)` before
#' running this function.
#'
#' In our paper, we report on the results of the following simulation study. For
#' each of the two models a) growth curve model and b) two-factor model, we
#' generate a sample size of `n` observations from a distribution `dist` with
#' reliability `rel`. We then fit the model using the maximum likelihood (`ML`)
#' estimator, the explicit RBM estimator (`eRBM`), the implicit RBM estimator
#' (`iRBM`), and the implicit RBM estimator with plug-in penalties (`iRBMp_x`).
#' We repeat this process `nsimu` times and report the results in terms of bias,
#' standard error, and computation time.
#'
#' @inheritParams gen-data
#' @inheritParams fit_sem
#' @param model The model to simulate. Either `growth` or `twofac`.
#' @param nsimu The number of replications to run.
#' @param whichsims A character vector with the estimators to use. Possible
#'   values are `ML`, `eRBM`, and `iRBM` for now
#'
#' @return A list with two elements. The first element is a data frame with the
#'   results of the simulation study. The second element is a list with the
#'   errors that occurred during the simulation study.
#' @export
sim_fun <- function(
    dist = "Normal",
    model = "twofac",
    rel = 0.8,
    n = 25,
    nsimu = 1,
    lavsim = FALSE,
    lavfun = "sem",
    whichsims = c("ML", "eRBM", "iRBM"),
    info_pen = "observed",
    info_bias = "observed",
    info_se = "observed"
  ) {
  dist <- match.arg(dist, c("Normal", "Kurtosis", "Non-normal"))
  model <- match.arg(model, c("growth", "twofac"))
  rel <- match.arg(as.character(rel), c("0.8", "0.5"))

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
  datasets <- replicate(
    nsimu,
    gen_data(n = n, rel = rel, dist = dist, lavsim = lavsim),
    simplify = FALSE
  )
  mod <- txt_mod(rel)
  true_vals <- truth(datasets[[1]])

  # Single run function --------------------------------------------------------
  single_sim <- function(j) {
    dat <- datasets[[j]]

    fitsemargs <- list(
      model = mod,
      data = dat,
      estimator = "ML",
      information = "expected",
      debug = FALSE,
      lavfun = lavfun,
      info_pen = info_pen,
      info_bias = info_bias,
      info_se = info_se,
      maxgrad = TRUE
    )

    nsimtypes <- length(whichsims)
    fit_list <- list()

    if ("ML" %in% whichsims) {
      fitsemargs$estimator.args <- list(rbm = "none", plugin_penalty = NULL)

      # Start from lavaan defaults
      fita <- do.call(fit_sem, fitsemargs)

      # Start from truth
      fitsemargs$start <- true_vals
      fitb <- do.call(fit_sem, fitsemargs)

      fitz <- list(fita, fitb)
      compare_loglik <- sapply(fitz, \(x) x$optim$objective)
      higher_loglik <- which.min(compare_loglik)

      fit_list$ML <- fitz[[higher_loglik]]
      fitsemargs$start <- fit_list$ML$coefficients  # use ML as starting values
    }
    if ("eRBM" %in% whichsims) {
      fitsemargs$estimator.args <- list(rbm = "explicit", plugin_penalty = NULL)
      fit_list$eRBM <- do.call(fit_sem, fitsemargs)
      fitsemargs$start <- fit_list$eRBM$coefficients
    }
    if ("iRBM" %in% whichsims) {
      fitsemargs$estimator.args <- list(rbm = "implicit", plugin_penalty = NULL)

      # Start from eBR estimates
      fitc <- do.call(fit_sem, fitsemargs)

      # Start from lavaan defaults
      fitsemargs$start <- NULL
      fita <- do.call(fit_sem, fitsemargs)

      # Start from truth
      fitsemargs$start <- true_vals
      fitb <- do.call(fit_sem, fitsemargs)

      # Start from ML estimates
      fitsemargs$start <- fit_list$ML$coefficients
      fitd <- do.call(fit_sem, fitsemargs)

      fitz <- list(fita, fitb, fitc, fitd)
      compare_loglik <- sapply(fitz, \(x) x$optim$objective)
      higher_loglik <- which.min(compare_loglik)
      fit_list$iRBM <- fitz[[higher_loglik]]
    }

    tibble::tibble(
      sim = j,
      dist = dist,
      model = model,
      rel = rel,
      n = n,
      info_pen = substr(info_pen, 1, 3),
      info_bias = substr(info_bias, 1, 3),
      info_se = substr(info_se, 1, 3),
      method = names(fit_list),
      est = lapply(fit_list, \(x) x$coefficients),
      se = lapply(fit_list, \(x) x$stderr),
      truth = rep(list(true_vals), nsimtypes),
      timing = sapply(fit_list, \(x) x$timing),
      converged = sapply(fit_list, \(x) x$converged),
      scaled_grad = lapply(fit_list, \(x) x$scaled_grad),
      max_loglik = sapply(fit_list, \(x) -1 * x$optim$objective),
      Sigma_OK = sapply(fit_list, \(x) {
        EV <- eigen(x$Sigma, only.values = TRUE)$values
        all(EV > 0)
      }),
      optim_message = sapply(fit_list, \(x) x$optim$message)
    )
  }

  # Run simulation -------------------------------------------------------------
  simu_res <- furrr::future_map(
    seq_len(nsimu),
    purrr::safely(single_sim, otherwise = NA),
    .progress = TRUE,
    .options = furrr::furrr_options(seed = TRUE)
  )

  # Clean up -------------------------------------------------------------------
  where_error <- which(sapply(simu_res, \(x) !is.null(x$error)))
  error_list <- lapply(simu_res[where_error], \(x) x$error)
  names(error_list) <- where_error

  # Return results -------------------------------------------------------------
  list(
    simu_res = do.call("rbind", lapply(simu_res, \(x) x$result)),
    errors = error_list
  )
}
