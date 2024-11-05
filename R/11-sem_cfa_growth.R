brcfa <- cfa <- function(
    model,
    data,
    estimator = c("iBRM", "iBRMp", "eBRM", "ML"),
    information = c("expected", "observed", "first.order"),
    ...
  ) {

  fit <- fit_sem(
    model = model,
    data = data,
    estimator = match.arg(estimator),
    information = match.arg(information),
    debug = FALSE,
    lavfun = "cfa",
    ...
  )

  # Prepare output
  x <- fit$coef
  lavargs <- list(...)
  lavargs$model <- model
  lavargs$data <- data
  lavargs$do.fit <- FALSE
  lavargs$information <- fit$information_se
  lavargs$start <- x
  fit0 <- do.call(get(fit$lavfun, envir = asNamespace("lavaan")), lavargs)

  # Change Model and implied slots
  fit0@Model <- lav_model_set_parameters(fit0@Model, x)
  fit0@implied <- lav_model_implied(fit0@Model)

  # Change ParTable and pta slots
  pt <- partable(fit0)
  pt$est[pt$free > 0] <- x
  pt$se <- 0
  pt$se[pt$free > 0] <- fit$stderr
  fit0@ParTable <- as.list(pt)
  fit0@pta$names <- names(pt)

  # Change Options slot
  fit0@Options$estimator <- "E-BRML"
  fit0@Options$estimator.args <- list(method = "eRBM")
  # fit0@Options$test <- "standard"
  fit0@Options$se <- "standard"

  # Change Fit slot
  fit0@Fit@x <- x
  fit0@Fit@se <- fit$stderr
  fit0@Fit@iterations <- fit$optim$iterations
  fit0@Fit@converged <- fit$optim$convergence == 0L

  # Change optim slot
  fit0@optim$x <- x
  # fit0@optim$dx <- 0
  fit0@optim$npar <- length(x)
  fit0@optim$fx <- fit0@Fit@fx
  fit0@optim$fx.group <- fit0@Fit@fx.group
  fit0@optim$iterations <- fit$optim$iterations
  fit0@optim$converged <- fit$optim$convergence == 0L

  # Change loglik slot
  fit0@loglik$estimator <- "E-BRML"

  # Change vcov slot
  fit0@vcov$se <- "standard"
  fit0@vcov$vcov <- fit$vcov

  # fit0@test <- fit_lav@test
  # fit0@baseline <- fit_lav@baseline

  fit0
}
