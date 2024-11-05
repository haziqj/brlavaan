sim_fun <- function(
    dist = "Normal",
    model = "twofac",
    rel = 0.8,
    n = 100,
    lavsim = FALSE,
    lavfun = "sem",
    meanstructure = FALSE,
    bounds = "none"
  ) {
  dist <- match.arg(dist, c("Normal", "Kurtosis", "Non-normal"))
  model <- match.arg(model, c("growth", "twofac"))
  rel <- match.arg(as.character(rel), c("0.8", "0.5"))

  # Generate all data ----------------------------------------------------------
  if (model == "growth") {
    gen_data <- gen_data_growth
    txt_mod <- txt_mod_growth
  } else {
    gen_data <- gen_data_twofac
    txt_mod <- txt_mod_twofac
  }
  datasets <- replicate(
    nsimu,
    gen_data(n = n, rel = rel, dist = dist, lavsim = lavsim, meanstructure = meanstructure),
    simplify = FALSE
  )
  mod <- txt_mod(rel)
  true_vals <- truth(datasets[[1]])

  simu_res <- future_map(
    seq_len(nsimu),
    safely(function(j) {
      dat <- datasets[[j]]

      fitsemargs <- list(
        model = mod,
        data = dat,
        lavfun = "cfa",
        meanstructure = meanstructure,
        bounds = "standard"
      )
      fitsemargs$method <- "ML"   ; fit_ML    <- do.call(fit_sem, fitsemargs)
      fitsemargs$method <- "eRBM" ; fit_eRBM  <- do.call(fit_sem, fitsemargs)
      fitsemargs$method <- "iRBM" ; fit_iRBM  <- do.call(fit_sem, fitsemargs)
      fitsemargs$method <- "iRBMp"; fit_iRBMp <- do.call(fit_sem, fitsemargs)

      data.frame(
        dist = dist,
        model = model,
        rel = rel,
        n = n,
        param = rep(names(true_vals), 4),
        est = c(fit_ML$coef, fit_eRBM$coef, fit_iRBM$coef, fit_iRBMp$coef),
        se = c(fit_ML$stderr, fit_eRBM$stderr, fit_iRBM$stderr, fit_iRBMp$stderr),
        truth = rep(true_vals, 4),
        method = rep(c("ML", "eRBM", "iRBM", "iRBMp"), each = length(true_vals)),
        time = rep(c(fit_ML$time, fit_eRBM$time, fit_iRBM$time, fit_iRBMp$time), each = length(true_vals))
      )
    }, otherwise = NA),
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
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
