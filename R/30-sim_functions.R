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
#' @param info_pen Should be `"observed"`
#' @param info_bias Should be `"observed"`
#' @param info_se Should be `"observed"`
#' @param bounds Bounded estimation?
#' @param keep_going If `TRUE`, the simulation will continue until the desired
#'   `nsimu` runs are obtained.
#' @param data_scale A scaling factor for the data. Default is `1`. Usefor for
#'   the growth curve simulations that this is set to `1/10`.
#' @param seeds A vector of seeds to use for the simulation study.
#'
#' @return A list with two elements. The first element is a data frame with the
#'   results of the simulation study. The second element is a list with the
#'   errors that occurred during the simulation study.
#' @export
sim_fun <- function(
    dist = "Normal",
    model = "twofac",
    rel = 0.8,
    n = 15,
    nsimu = 1,
    lavsim = FALSE,
    whichsims = c("ML", "eRBM", "iRBM"),
    bounds = "none",
    info_pen = "observed",
    info_bias = "observed",
    info_se = "observed",
    keep_going = FALSE,
    data_scale = 1,
    seeds = NULL
  ) {
  dist <- match.arg(dist, c("Normal", "Kurtosis", "Non-normal"))
  model <- match.arg(model, c("growth", "twofac"))
  rel <- match.arg(as.character(rel), c("0.8", "0.5"))

  # Check seeds
  if (!is.null(seeds)) {
    if (length(seeds) != nsimu) {
      cli::cli_abort("Length of seeds must be equal to nsimu")
    }
    if (isTRUE(keep_going)) {
      cli::cli_abort("Cannot use keep_going and seeds at the same time")
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
      info_pen = info_pen,
      info_bias = info_bias,
      info_se = info_se,
      debug = FALSE,
      lavfun = lavfun,
      maxgrad = TRUE,
      nearPD = FALSE,

      bounds = bounds,
      start = true_vals
    )

    nsimtypes <- length(whichsims)
    fit_list <- list()

    if ("ML" %in% whichsims) {
      fit_list$ML <- try(do.call(fit_sem, fitsemargs))
    }
    if ("eRBM" %in% whichsims) {
      fitsemargs$rbm <- "explicit"
      fit_list$eRBM <- try(do.call(fit_sem, fitsemargs))
    }
    if ("iRBM" %in% whichsims) {
      fitsemargs$rbm <- "implicit"
      fit_list$iRBM <- do.call(fit_sem, fitsemargs)
    }

    # converged <- sapply(fit_list, \(x) x$converged)
    # if (any(!converged)) return(NULL)

    tibble::tibble(
      seed = seeds[j],
      sim = j,
      dist = dist,
      model = model,
      rel = rel,
      n = n,
      info_pen = substr(info_pen, 1, 3),
      info_bias = substr(info_bias, 1, 3),
      info_se = substr(info_se, 1, 3),
      method = names(fit_list),
      est = lapply(fit_list, purrr::possibly(\(x) x$coefficients, NA)),
      se = lapply(fit_list, purrr::possibly(\(x) x$stderr, NA)),
      truth = rep(list(true_vals), nsimtypes),
      timing = sapply(fit_list, purrr::possibly(\(x) x$timing, NA)),
      converged = sapply(fit_list, purrr::possibly(\(x) x$converged, NA)),
      scaled_grad = lapply(fit_list, purrr::possibly(\(x) x$scaled_grad, NA)),
      max_loglik = sapply(fit_list, purrr::possibly(\(x) -1 * x$optim$objective, NA)),
      Sigma_OK = sapply(fit_list, purrr::possibly(\(x) {
        EV <- eigen(x$Sigma, only.values = TRUE)$values
        all(EV > 0)
      }, NA)),
      optim_message = sapply(fit_list, purrr::possibly(\(x) x$optim$message, NA))
    )
  }

  # Run simulation -------------------------------------------------------------
  if (!isFALSE(keep_going)) {
    simu_res <- list()
    i <- 1
    how_many_fine <- 0
    while (how_many_fine < nsimu) {
      batch_size <- nsimu - how_many_fine
      cat("\n")
      cli::cli_alert_info("Batch size: {batch_size} \n")
      new_res <- furrr::future_map(
        .x = (i - 1) + seq_len(batch_size),
        .f = purrr::safely(single_sim, otherwise = NULL),
        .progress = TRUE,
        .options = furrr::furrr_options(seed = TRUE)
      )

      # Convergence table
      conv_count <- sum(sapply(new_res, \(x) {
        y <- x$result$converged  # vector of length 3 in order of ML, eRBM, iRBM
        if (is.null(y)) return(FALSE)
        else {
          if (is.numeric(keep_going)) {
            out <- all(y[keep_going])
          } else {
            out <- all(y)
          }
          if (is.na(out)) out <- FALSE
          return(out)
        }
      }))
      how_many_fine <- how_many_fine + conv_count

      simu_res <- append(simu_res, new_res)
      i <- i + batch_size
    }
  } else {
    simu_res <- furrr::future_map(
      seq_len(nsimu),
      purrr::safely(single_sim, otherwise = NA),
      .progress = TRUE,
      .options = furrr::furrr_options(seed = TRUE)
    )
  }

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
