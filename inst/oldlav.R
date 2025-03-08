library(brlavaan)
i <- 5
check <- TRUE
while (check) {
  i <- i + 1
  set.seed(i)
  dat <- gen_data_growth(n = 15, rel = 0.5, dist = "Normal", scale = 1/10)
  mod <- txt_mod_growth(0.5)
  tru <- truth(dat)

  fit_lav    <- growth(mod, dat, start = tru)
  fit_ML     <- fit_sem(mod, dat, lavfun = "growth", start = tru, rbm = "none")
  check <- tinytest::expect_equal(
    as.numeric(coef(fit_lav)),
    as.numeric(coef(fit_ML)),
    tolerance = 1e-3
  )
  print(i)
}

print(check)
plot(coef(fit_lav), coef(fit_ML))
abline(0, 1)

# library(brlavaan)
# rel <- 0.5
# mod <- txt_mod_growth(rel)
# dat <- gen_data_growth(n = 15, rel, dist = "Normal")
# cov(dat)
#
#
# remove.packages("lavaan")
# devtools::install_version("lavaan", "0.6-9")  # 2021-06-27
# library(lavaan)
# tinytest::expect_equal(packageVersion("lavaan"), package_version("0.6-9"))
# fit_old <- growth(mod, dat, meanstructure = TRUE, bounds = "standard")
# coef(fit_old)
#
# remove.packages("lavaan")
# devtools::install_cran("lavaan")
# library(lavaan)
# tinytest::expect_equal(packageVersion("lavaan"), package_version("0.6-19"))
# fit_new <- growth(mod, dat, meanstructure = TRUE, bounds = "standard")
# coef(fit_new)
#
#
