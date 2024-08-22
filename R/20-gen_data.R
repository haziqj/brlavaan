# Growth curve models (model a) ------------------------------------------------
txt_mod_growth <- function(rel) {
  if (rel == "0.8") {
    # Mean reliability 0.80 ----------------------------------------------------
    mod <- "
    i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4 + 1*y5 + 1*y6 + 1*y7 + 1*y8 + 1*y9 + 1*y10
    s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4 + 4*y5 + 5*y6 + 6*y7 + 7*y8 + 8*y9 + 9*y10

    i ~~ 550 * i
    s ~~ 100 * s
    i ~~ 40 * s
    y1 ~~ 500 * y1
    y2 ~~ 500 * y2
    y3 ~~ 500 * y3
    y4 ~~ 500 * y4
    y5 ~~ 500 * y5
    y6 ~~ 500 * y6
    y7 ~~ 500 * y7
    y8 ~~ 500 * y8
    y9 ~~ 500 * y9
    y10 ~~ 500 * y10
    "
  }
  if (rel == "0.5") {
    # Mean reliability 0.50 ----------------------------------------------------
    mod <- "
    i =~ 1*y1 + 1*y2 + 1*y3 + 1*y4 + 1*y5 + 1*y6 + 1*y7 + 1*y8 + 1*y9 + 1*y10
    s =~ 0*y1 + 1*y2 + 2*y3 + 3*y4 + 4*y5 + 5*y6 + 6*y7 + 7*y8 + 8*y9 + 9*y10

    i ~~ 275 * i
    s ~~ 50 * s
    i ~~ 20 * s
    y1 ~~ 1300 * y1
    y2 ~~ 1300 * y2
    y3 ~~ 1300 * y3
    y4 ~~ 1300 * y4
    y5 ~~ 1300 * y5
    y6 ~~ 1300 * y6
    y7 ~~ 1300 * y7
    y8 ~~ 1300 * y8
    y9 ~~ 1300 * y9
    y10 ~~ 1300 * y10
    "
  }

  mod
}

truth_growth <- function(rel) {
  if (rel == "0.8") {
    truth <- c(rep(500, 10), 550, 100, 40)
  }
  if (rel == "0.5") {
    truth <- c(rep(1300, 10), 275, 50, 20)
  }

  truth
}

gen_data_growth <- function(n = 100, rel = 0.8) {
  mod <- txt_mod_growth(rel = as.character(rel))
  # truth <- truth_growth(rel)
  dat <- simulateData(model = mod, sample.nobs = n)
  # dat$truth <- truth
  dat
}

# Two-factor SEM models (model b) ----------------------------------------------
txt_mod_twofac <- function(rel) {
  if (rel == "0.8") {
    # Reliabilities 0.80 -------------------------------------------------------
    mod <- "
    eta1 =~ 1*y1 + 0.7*y2 + 0.6*y3
    eta2 =~ 1*y4 + 0.7*y5 + 0.6*y6
    eta2 ~ 0.25*eta1

    eta1 ~~ 1*eta1
    eta2 ~~ 1*eta2
    y1 ~~ 0.25*y1
    y2 ~~ 0.09*y2
    y3 ~~ 0.1225*y3
    y4 ~~ 0.25*y4
    y5 ~~ 0.09*y5
    y6 ~~ 0.1225*y6
    "
  }
  if (rel == "0.5") {
    # Reliabilities 0.50 -------------------------------------------------------
    mod <- "
    eta1 =~ 1*y1 + 0.7*y2 + 0.6*y3
    eta2 =~ 1*y4 + 0.7*y5 + 0.6*y6
    eta2 ~ 0.25*eta1

    eta1 ~~ 1*eta1
    eta2 ~~ 1*eta2
    y1 ~~ 1*y1
    y2 ~~ 0.49*y2
    y3 ~~ 0.36*y3
    y4 ~~ 1*y4
    y5 ~~ 0.49*y5
    y6 ~~ 0.36*y6
    "
  }

  mod
}

truth_twofac <- function(rel) {
  if (rel == "0.8") {
    truth <- c(0.7, 0.6, 0.7, 0.6, 0.25, rep(c(0.25, 0.09, 0.1225), 2), 1, 1)
  }
  if (rel == "0.5") {
    truth <- c(0.7, 0.6, 0.7, 0.6, 0.25, rep(c(1, 0.49, 0.36), 2), 1, 1)
  }

  truth
}

gen_data_twofac <- function(n = 100, rel = 0.8) {
  mod <- txt_mod_twofac(rel = as.character(rel))
  # truth <- truth_twofac(rel)
  dat <- simulateData(model = mod, sample.nobs = n)
  # dat$truth <- truth
  dat
}
