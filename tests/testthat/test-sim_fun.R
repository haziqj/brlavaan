test_that("Two-factor simulation function works", {
  expect_no_error({
    future::plan(future::sequential)
    res <- sim_fun(
      dist = "Normal",
      model = "twofac",
      rel = 0.8,
      n = 100,
      nsimu = 1,
      lavsim = TRUE,
      lavfun = "sem"
    )
  })
})

test_that("Growth simulation function works", {
  expect_no_error({
    future::plan(future::sequential)
    res <- sim_fun(
      dist = "Kurtosis",
      model = "growth",
      rel = 0.5,
      n = 50,
      nsimu = 1,
      lavsim = FALSE,
      lavfun = "growth"
    )
  })
})
