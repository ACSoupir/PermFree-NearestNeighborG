test_that("univariate exact moments match exhaustive enumeration", {
  set.seed(42)
  N <- 10L
  n <- 4L
  x <- stats::runif(N)
  y <- stats::runif(N)
  radii <- seq(0, 1.2, by = 0.15)
  D <- as.matrix(stats::dist(cbind(x, y)))
  M <- enum_uni(D, radii, n)

  cf <- csr_moments_r(x, y, radii, n_i = n)
  ccpp <- permfreeG:::csr_uni_moments_cpp(x, y, radii, n)

  expect_equal(cf$mean, pop_mean(M), tolerance = 1e-10)
  expect_equal(cf$var, pop_var(M), tolerance = 1e-10)
  expect_equal(ccpp$mean, pop_mean(M), tolerance = 1e-10)
  expect_equal(ccpp$var, pop_var(M), tolerance = 1e-10)
})

test_that("public API reproduces the exact mean and variance", {
  set.seed(42)
  N <- 10L
  n <- 4L
  d <- data.frame(x = stats::runif(N), y = stats::runif(N),
                  m = factor(c(rep("A", n), rep("B", N - n))))
  radii <- seq(0, 1.2, by = 0.15)
  res <- as.data.frame(exact_gest(d, radii = radii, mark_col = "m",
                                  marks_i = "A"))
  D <- as.matrix(stats::dist(cbind(d$x, d$y)))
  M <- enum_uni(D, radii, n)
  expect_equal(res$csr_mean, pop_mean(M), tolerance = 1e-10)
  expect_equal(res$csr_var, pop_var(M), tolerance = 1e-10)
  expect_true(all(is.na(res$z[res$csr_var == 0])))
})

test_that("independence approximation changes only the variance", {
  set.seed(7)
  N <- 20L
  n <- 6L
  x <- stats::runif(N)
  y <- stats::runif(N)
  radii <- seq(0, 1.5, by = 0.25)
  full <- permfreeG:::csr_uni_moments_cpp(x, y, radii, n)
  indep <- permfreeG:::csr_uni_moments_cpp(x, y, radii, n, covariance = FALSE)
  expect_equal(full$mean, indep$mean, tolerance = 1e-12)
  expect_true(all(full$var >= 0))
  expect_true(all(indep$var >= 0))
})

test_that("exact mean agrees with a permutation run", {
  skip_on_cran()
  set.seed(1)
  N <- 20L
  n <- 6L
  d <- data.frame(x = stats::runif(N), y = stats::runif(N),
                  m = factor(c(rep("A", n), rep("B", N - n))))
  radii <- seq(0, 1.5, by = 0.25)
  ex <- as.data.frame(exact_gest(d, radii = radii, mark_col = "m",
                                 marks_i = "A"))
  pm <- permute_gest(d, radii = radii, mark_col = "m", marks_i = "A",
                     nsim = 400L)
  expect_equal(ex$csr_mean, pm$mean, tolerance = 0.05)
})
