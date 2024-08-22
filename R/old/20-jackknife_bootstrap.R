# FIXME: Note that the below uses semList(). Maybe in the future want to do
# lavaanList() with optional cmd argument for sem/cfa/efa etc.

# Jackknife bias correction ----------------------------------------------------
rb_jack <- function(fit, data) {

  start_time <- Sys.time()

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

  end_time <- Sys.time()

  # Compute bias-corrected estimates (Equation 10)
  out <- nrow(data) * coef(fit) -
    (nrow(data)- 1 ) * rowMeans(coef(jack.estimates), na.rm = TRUE)

  attr(out, "meta") <- jack.estimates@meta
  attr(out, "timing") <- end_time - start_time
  out
}

# Bootstrap bias correction ----------------------------------------------------
rb_boot <- function(fit, data, bootsamp = 500) {

  start_time <- Sys.time()

  resample.bootstrap <- function(.data){
    df.list <- vector(mode = "list", length = bootsamp)
    for(i in seq_len(bootsamp)) {
      rowz <- sample.int(nrow(.data), size = nrow(.data), replace = TRUE)
      df.list[[i]] <- .data[rowz, ]
    }
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

  end_time <- Sys.time()

  # Compute bias-corrected estimates (Equation 12)
  out <- 2 * coef(fit) - rowMeans(coef(boot.estimates), na.rm = TRUE)

  attr(out, "meta") <- boot.estimates@meta
  attr(out, "timing") <- end_time - start_time
  out
}
