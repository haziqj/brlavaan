# Jackknife bias correction ----------------------------------------------------
rb_jack <- function(fit, data) {
  resample.jackknife <- function(.data){
    df.list <- vector(mode = "list", length = nrow(.data))
    for(i in 1:nrow(.data)) { df.list[[i]] <- .data[-i, ] }
    return(df.list)
  }

  # Create a list of all n jackknife resamples
  jack.resamples <- resample.jackknife(data)

  # Fit the model to all dataframes in the list
  jack.estimates <- semList(
    model = rlang::call_args(fit@call)$model,
    dataList = jack.resamples,
    store.slots = "partable"
  )

  # Compute bias-corrected estimates (Equation 10)
  nrow(data) * coef(fit) - (nrow(data)- 1 ) * rowMeans(coef(jack.estimates))
}

# Bootstrap bias correction ----------------------------------------------------
rb_boot <- function(fit, data) {
  resample.bootstrap <- function(.data){
    df.list <- vector(mode = "list", length = 500)
    for(i in 1:500) { df.list[[i]] <- data[sample.int(nrow(.data),
                                                      size = nrow(.data),
                                                      replace = TRUE), ] }
    return(df.list)
  }

  # Create a list of all 500 Bootstrap resamples
  boot.resamples <- resample.bootstrap(data)

  # Fit the model to all dataframes in the list
  boot.estimates <- semList(
    model = rlang::call_args(fit@call)$model,
    dataList = boot.resamples,
    store.slots = "partable"
  )

  # Compute bias-corrected estimates (Equation 12)
  2 * coef(fit) - rowMeans(coef(boot.estimates))
}
