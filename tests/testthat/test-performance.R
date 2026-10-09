test_that("exact engine handles N=250 x 151 radii quickly", {
  skip_on_cran()
  set.seed(1)
  N <- 250L
  x <- stats::runif(N)
  y <- stats::runif(N)
  radii <- seq(0, 1.5, length.out = 151L)
  expect_under({
    m <- permfreeG:::csr_uni_moments_cpp(x, y, radii, 50L)
    expect_length(m$var, 151L)
    expect_true(all(is.finite(m$mean)))
  }, secs = 30)
})

test_that("bivariate exact engine handles N=250 x 101 radii quickly", {
  skip_on_cran()
  set.seed(1)
  N <- 250L
  x <- stats::runif(N)
  y <- stats::runif(N)
  radii <- seq(0, 1.5, length.out = 101L)
  expect_under({
    m <- permfreeG:::csr_cross_moments_cpp(x, y, radii, 40L, 35L)
    expect_length(m$var, 101L)
    expect_true(all(is.finite(m$mean)))
  }, secs = 45)
})

test_that("approximation is faster than the exact variance", {
  skip_on_cran()
  set.seed(1)
  N <- 250L
  x <- stats::runif(N)
  y <- stats::runif(N)
  radii <- seq(0, 1.5, length.out = 101L)
  t_exact <- system.time(permfreeG:::csr_uni_moments_cpp(x, y, radii, 50L))[["elapsed"]]
  t_approx <- system.time(
    permfreeG:::csr_uni_moments_cpp(x, y, radii, 50L, covariance = FALSE))[["elapsed"]]
  expect_true(t_approx <= t_exact + 0.05)
})
