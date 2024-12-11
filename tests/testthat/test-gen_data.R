## ----- Growth curve model ----------------------------------------------------
test_that("Normal distribution", {
  set.seed(123)
  dat1 <- gen_data_growth(n = 5, rel = 0.8, dist = "Normal", lavsim = FALSE)
  dat2 <- gen_data_growth(n = 5, rel = 0.5, dist = "Normal", lavsim = FALSE)

  # expect_snapshot(dat1)
  # expect_snapshot(dat2)
  expect_equal(truth(dat1), truth_growth(0.8))
})

test_that("lavaan's normal distribution", {
  set.seed(123)
  dat1 <- gen_data_growth(n = 5, rel = 0.8, lavsim = TRUE)
  dat2 <- gen_data_growth(n = 5, rel = 0.5, lavsim = FALSE)

  # expect_snapshot(dat1)
  # expect_snapshot(dat2)
  expect_equal(truth(dat2), truth_growth(0.5))
})

test_that("Kurtosis distribution", {
  set.seed(123)
  dat1 <- gen_data_growth(n = 5, rel = 0.8, dist = "Kurtosis", lavsim = FALSE)
  dat2 <- gen_data_growth(n = 5, rel = 0.5, dist = "Kurtosis", lavsim = FALSE)

  # expect_snapshot(dat1)
  # expect_snapshot(dat2)
  expect_equal(truth(dat1), truth_growth(0.8))
})

test_that("Non-normal distribution", {
  set.seed(123)
  dat1 <- gen_data_growth(n = 5, rel = 0.8, dist = "Non-normal", lavsim = FALSE)
  dat2 <- gen_data_growth(n = 5, rel = 0.5, dist = "Non-normal", lavsim = FALSE)

  # expect_snapshot(dat1)
  # expect_snapshot(dat2)
  expect_equal(truth(dat2), truth_growth(0.5))
})

## ---- Two-factor model -------------------------------------------------------
test_that("Normal distribution", {
  set.seed(123)
  dat1 <- gen_data_twofac(n = 5, rel = 0.8, dist = "Normal", lavsim = FALSE)
  dat2 <- gen_data_twofac(n = 5, rel = 0.5, dist = "Normal", lavsim = FALSE)

  # expect_snapshot(dat1)
  # expect_snapshot(dat2)
  expect_equal(truth(dat1), truth_twofac(0.8))
})

test_that("lavaan's normal distribution", {
  set.seed(123)
  dat1 <- gen_data_twofac(n = 5, rel = 0.8, lavsim = TRUE)
  dat2 <- gen_data_twofac(n = 5, rel = 0.5, lavsim = FALSE)

  # expect_snapshot(dat1)
  # expect_snapshot(dat2)
  expect_equal(truth(dat2), truth_twofac(0.5))
})

test_that("Kurtosis distribution", {
  set.seed(123)
  dat1 <- gen_data_twofac(n = 5, rel = 0.8, dist = "Kurtosis", lavsim = FALSE)
  dat2 <- gen_data_twofac(n = 5, rel = 0.5, dist = "Kurtosis", lavsim = FALSE)

  # expect_snapshot(dat1)
  # expect_snapshot(dat2)
  expect_equal(truth(dat1), truth_twofac(0.8))
})

test_that("Non-normal distribution", {
  set.seed(123)
  dat1 <- gen_data_twofac(n = 5, rel = 0.8, dist = "Non-normal", lavsim = FALSE)
  dat2 <- gen_data_twofac(n = 5, rel = 0.5, dist = "Non-normal", lavsim = FALSE)

  # expect_snapshot(dat1)
  # expect_snapshot(dat2)
  expect_equal(truth(dat2), truth_twofac(0.5))
})
