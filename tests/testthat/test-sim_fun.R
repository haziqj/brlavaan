test_that("Two-factor simulation function works", {
  future::plan(future::sequential)
  res <- sim_fun(
    dist = "Normal",
    model = "twofac",
    rel = 0.8,
    n = 100,
    nsimu = 1,
    lavsim = TRUE
  )

  expect_type(res, "list")
  expect_s3_class(res$simu_res, "data.frame")
})

test_that("Growth simulation function works", {
  future::plan(future::sequential)
  res <- sim_fun(
    dist = "Kurtosis",
    model = "growth",
    rel = 0.5,
    n = 50,
    nsimu = 1,
    lavsim = FALSE
  )

  expect_type(res, "list")
  expect_s3_class(res$simu_res, "data.frame")
})

test_that("Seeds work", {
  future::plan(future::sequential)

  expect_error(sim_fun(nsimu = 10, seeds = 1234))

  res1 <- sim_fun(
    dist = "Normal",
    model = "twofac",
    rel = 0.8,
    n = 100,
    nsimu = 1,
    lavsim = TRUE,
    whichsims = "ML",
    seeds = 1234
  )
  res2 <- sim_fun(
    dist = "Normal",
    model = "twofac",
    rel = 0.8,
    n = 100,
    nsimu = 1,
    lavsim = TRUE,
    whichsims = "ML",
    seeds = 1234
  )

  expect_equal(res1$simu_res$est, res2$simu_res$est)

})
