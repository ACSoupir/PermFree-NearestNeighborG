test_that("bivariate exact moments match exhaustive enumeration", {
  set.seed(11)
  N <- 10L
  ni <- 4L
  nj <- 3L
  x <- stats::runif(N)
  y <- stats::runif(N)
  radii <- seq(0, 1.2, by = 0.15)
  D <- as.matrix(stats::dist(cbind(x, y)))
  M <- enum_cross(D, radii, ni, nj)

  cf <- csr_moments_r(x, y, radii, n_i = ni, n_j = nj)
  ccpp <- permfreeG:::csr_cross_moments_cpp(x, y, radii, ni, nj)

  expect_equal(cf$mean, pop_mean(M), tolerance = 1e-10)
  expect_equal(cf$var, pop_var(M), tolerance = 1e-9)
  expect_equal(ccpp$mean, pop_mean(M), tolerance = 1e-10)
  expect_equal(ccpp$var, pop_var(M), tolerance = 1e-9)
})

test_that("cross mean depends on n_j but not on n_i", {
  set.seed(5)
  N <- 30L
  x <- stats::runif(N)
  y <- stats::runif(N)
  radii <- seq(0, 1.5, by = 0.3)
  a <- permfreeG:::csr_cross_moments_cpp(x, y, radii, 5L, 8L)
  b <- permfreeG:::csr_cross_moments_cpp(x, y, radii, 9L, 8L)
  cc <- permfreeG:::csr_cross_moments_cpp(x, y, radii, 5L, 12L)
  expect_equal(a$mean, b$mean, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(a$mean, cc$mean)))
})

test_that("public bivariate API matches exhaustive enumeration", {
  set.seed(11)
  N <- 10L
  ni <- 4L
  nj <- 3L
  d <- data.frame(x = stats::runif(N), y = stats::runif(N),
                  m = factor(c(rep("A", ni), rep("B", nj),
                               rep("C", N - ni - nj))))
  radii <- seq(0, 1.2, by = 0.15)
  res <- as.data.frame(exact_gcross(d, radii = radii, mark_col = "m",
                                    marks_i = "A", marks_j = "B"))
  D <- as.matrix(stats::dist(cbind(d$x, d$y)))
  M <- enum_cross(D, radii, ni, nj)
  expect_equal(res$csr_mean, pop_mean(M), tolerance = 1e-10)
  expect_equal(res$csr_var, pop_var(M), tolerance = 1e-9)
})

test_that("i == j reduces to the univariate Gest statistic", {
  radii <- seq(0, 1.2, by = 0.3)
  gij <- exact_gcross(sim_nnG, radii = radii, mark_col = "m",
                      marks_i = "A", marks_j = "A")
  gi <- exact_gest(sim_nnG, radii = radii, mark_col = "m", marks_i = "A")
  expect_identical(gij$statistic, "Gest")
  a <- as.data.frame(gij)
  b <- as.data.frame(gi)
  expect_equal(a$csr_mean, b$csr_mean, tolerance = 1e-12)
  expect_equal(a$obs, b$obs, tolerance = 1e-12)
})

test_that("bivariate mean agrees with a permutation run", {
  skip_on_cran()
  set.seed(2)
  N <- 20L
  ni <- 5L
  nj <- 4L
  d <- data.frame(x = stats::runif(N), y = stats::runif(N),
                  m = factor(c(rep("A", ni), rep("B", nj),
                               rep("C", N - ni - nj))))
  radii <- seq(0, 1.5, by = 0.25)
  ex <- as.data.frame(exact_gcross(d, radii = radii, mark_col = "m",
                                   marks_i = "A", marks_j = "B"))
  pm <- permute_gcross(d, radii = radii, mark_col = "m",
                       marks_i = "A", marks_j = "B", nsim = 200L)
  expect_equal(ex$csr_mean, pm$mean, tolerance = 0.08)
})
