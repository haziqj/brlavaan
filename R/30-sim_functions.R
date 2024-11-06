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
#' estimator, the explicit BR ML estimator (`eBRM`), the implicit BR ML
#' estimator (`iBRM`), and the implicit BR ML estimator with a plug-in penalty
#' (`iBRMp`). We repeat this process `nsimu` times and report the results in
#' terms of bias, standard error, and computation time.
#'
#' @inheritParams gen-data
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
        lavfun = lavfun
      )
      fitsemargs$estimator <- "ML"   ; fit_ML    <- do.call(fit_sem, fitsemargs)
      fitsemargs$estimator <- "eBRM" ; fit_eBRM  <- do.call(fit_sem, fitsemargs)
      fitsemargs$estimator <- "iBRM" ; fit_iBRM  <- do.call(fit_sem, fitsemargs)
      fitsemargs$estimator <- "iBRMp"; fit_iBRMp <- do.call(fit_sem, fitsemargs)

      tibble::tibble(
        simu = j,
        dist = dist,
        model = model,
        rel = rel,
        n = n,
        param = rep(names(true_vals), 4),
        est = c(
          fit_ML$coefficients,
          fit_eBRM$coefficients,
          fit_iBRM$coefficients,
          fit_iBRMp$coefficients
        ),
        se = c(
          fit_ML$stderr,
          fit_eBRM$stderr,
          fit_iBRM$stderr,
          fit_iBRMp$stderr
        ),
        truth = rep(true_vals, 4),
        estimator = rep(
          c("ML", "eBRM", "iBRM", "iBRMp"),
          each = length(true_vals)
        ),
        timing = rep(
          c(fit_ML$timing, fit_eBRM$timing, fit_iBRM$timing, fit_iBRMp$timing),
          each = length(true_vals)
        ),
        converged = rep(
          c(fit_ML$converged, fit_eBRM$converged, fit_iBRM$converged,
            fit_iBRMp$converged),
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
