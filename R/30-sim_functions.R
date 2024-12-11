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
#' estimator, the explicit RBM estimator (`eRBM`), the implicit RBM
#' estimator (`iRBM`), and the implicit RBM estimator with plug-in penalties
#' (`iRBMp_x`). We repeat this process `nsimu` times and report the results in
#' terms of bias, standard error, and computation time.
#'
#' @inheritParams gen-data
#' @inheritParams fit_sem
#' @param model The model to simulate. Either `growth` or `twofac`.
#' @param nsimu The number of replications to run.
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
    nsimu = 4,
    lavsim = FALSE,
    lavfun = "sem"
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

  simu_res <- furrr::future_map(
    seq_len(nsimu),
    purrr::safely(function(j) {
      dat <- datasets[[j]]

      fitsemargs <- list(
        model = mod,
        data = dat,
        estimator = "ML",
        information = "expected",
        debug = FALSE,
        lavfun = lavfun
      )

      # ML
      fitsemargs$estimator.args <- list(rbm = "none", plugin_penalty = NULL)
      fit_ML    <- do.call(fit_sem, fitsemargs)

      # eRBM
      fitsemargs$estimator.args <- list(rbm = "explicit", plugin_penalty = NULL)
      fit_eRBM  <- do.call(fit_sem, fitsemargs)

      # iRBM
      fitsemargs$estimator.args <- list(rbm = "implicit", plugin_penalty = NULL)
      fit_iRBM  <- do.call(fit_sem, fitsemargs)

      # iRBMp ridge
      fitsemargs$estimator.args <- list(rbm = "implicit", plugin_penalty = pen_ridge)
      fit_iRBMpridge <- do.call(fit_sem, fitsemargs)

      # iRBMp ridge bounded
      fitsemargs$estimator.args <- list(rbm = "implicit", plugin_penalty = pen_ridge_bound)
      fit_iRBMpridgeb <- do.call(fit_sem, fitsemargs)

      # iRBMp Huber
      fitsemargs$estimator.args <- list(rbm = "implicit", plugin_penalty = pen_huber)
      fit_iRBMphuber <- do.call(fit_sem, fitsemargs)

      tibble::tibble(
        simu = j,
        dist = dist,
        model = model,
        rel = rel,
        n = n,
        param = rep(names(true_vals), 6),
        est = c(
          fit_ML$coefficients,
          fit_eRBM$coefficients,
          fit_iRBM$coefficients,
          fit_iRBMpridge$coefficients,
          fit_iRBMpridgeb$coefficients,
          fit_iRBMphuber$coefficients
        ),
        se = c(
          fit_ML$stderr,
          fit_eRBM$stderr,
          fit_iRBM$stderr,
          fit_iRBMpridge$stderr,
          fit_iRBMpridgeb$stderr,
          fit_iRBMphuber$stderr
        ),
        truth = rep(true_vals, 6),
        method = rep(
          c("ML", "eRBM", "iRBM", "iRBMp_ridge", "iRBMp_ridgebound", "iRBMp_huber"),
          each = length(true_vals)
        ),
        timing = rep(
          c(fit_ML$timing, fit_eRBM$timing, fit_iRBM$timing, fit_iRBMpridge$timing,
            fit_iRBMpridgeb$timing, fit_iRBMphuber$timing),
          each = length(true_vals)
        ),
        converged = rep(
          c(fit_ML$converged, fit_eRBM$converged, fit_iRBM$converged, fit_iRBMpridge$converged, fit_iRBMpridgeb$converged, fit_iRBMphuber$converged),
          each = length(true_vals)
        )
      )
    }, otherwise = NA),
    .progress = TRUE,
    .options = furrr::furrr_options(seed = TRUE)
  )
  simu_res

  # Clean up -------------------------------------------------------------------
  where_error <- which(sapply(simu_res, \(x) !is.null(x$error)))
  error_list <- lapply(simu_res[where_error], \(x) x$error)
  names(error_list) <- where_error

  list(
    simu_res = do.call("rbind", lapply(simu_res, \(x) x$result)),
    errors = error_list
  )
}
