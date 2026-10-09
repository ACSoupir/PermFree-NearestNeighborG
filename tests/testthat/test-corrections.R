test_that("border correction matches subset enumeration", {
  set.seed(3)
  N <- 10L
  n <- 4L
  d <- data.frame(x = stats::runif(N), y = stats::runif(N))
  d$m <- factor(c(rep("A", n), rep("B", N - n)))
  pp <- spatstat.geom::ppp(d$x, d$y,
                           window = spatstat.geom::owin(c(0, 1), c(0, 1)),
                           marks = d$m)
  radii <- seq(0, 0.3, by = 0.05)
  bd <- spatstat.geom::bdist.points(pp)
  D <- as.matrix(stats::dist(cbind(d$x, d$y)))
  ci <- utils::combn(N, n)
  est <- matrix(NA_real_, ncol(ci), length(radii))
  for (c in seq_len(ncol(ci))) {
    S <- ci[, c]
    nnd <- apply(D[S, S], 1, function(z) min(z[z > 0]))
    for (k in seq_along(radii)) {
      den <- sum(bd[S] > radii[k])
      est[c, k] <- if (den > 0) {
        sum(nnd <= radii[k] & bd[S] > radii[k]) / den
      } else NA_real_
    }
  }
  ok <- colSums(!is.na(est)) == ncol(ci)
  expect_true(sum(ok) >= 4L)

  gr <- as.data.frame(exact_gest(pp, radii = radii, marks_i = "A",
                                 correction = "rs"))
  em <- colMeans(est[, ok, drop = FALSE])
  ev <- apply(est[, ok, drop = FALSE], 2L,
              function(z) mean((z - mean(z))^2))
  # the reduced-sample mean is a ratio of expectations (first-order exact)
  expect_true(all(abs(gr$csr_mean[ok] - em) <= pmax(0.05 * abs(em), 1e-8)))
  # delta-method variance is a first-order approximation
  expect_true(all(abs(gr$csr_var[ok] - ev) <= pmax(0.35 * abs(ev), 1e-6)))
})

test_that("han and km mean-field curves track the permutation mean", {
  skip_on_cran()
  set.seed(99)
  N <- 40L
  n <- 15L
  d <- data.frame(x = stats::runif(N), y = stats::runif(N))
  d$m <- factor(c(rep("A", n), rep("B", N - n)))
  pp <- spatstat.geom::ppp(d$x, d$y,
                           window = spatstat.geom::owin(c(0, 1), c(0, 1)),
                           marks = d$m)
  radii <- seq(0, 0.4, by = 0.05)
  for (cor in c("han", "km")) {
    ex <- as.data.frame(exact_gest(pp, radii = radii, marks_i = "A",
                                   correction = cor))
    pm <- permute_gest(pp, radii = radii, marks_i = "A",
                       correction = cor, nsim = 200L)
    fin <- is.finite(ex$csr_mean) & is.finite(pm$mean)
    expect_true(any(fin))
    expect_equal(ex$csr_mean[fin], pm$mean[fin], tolerance = 0.12)
    expect_true(all(is.na(ex$csr_var)))
  }
})

test_that("han and km return well-formed curves", {
  radii <- seq(0, 1.5, by = 0.25)
  for (cor in c("han", "km")) {
    res <- exact_gest(sim_nnG, radii = radii, mark_col = "m", marks_i = "A",
                      correction = cor)
    df <- as.data.frame(res)
    expect_true(all(is.finite(df$csr_mean) | is.na(df$csr_mean)))
    fin <- df$csr_mean[is.finite(df$csr_mean)]
    expect_true(all(fin >= 0 & fin <= 1))
    expect_true(all(is.na(df$csr_var)))
    expect_false(res$exact)
  }
})

test_that("cross han/km fail with an informative message", {
  radii <- seq(0, 1.5, by = 0.25)
  expect_error(exact_gcross(sim_nnG, radii = radii, mark_col = "m",
                            marks_i = "A", marks_j = "B",
                            correction = "han"),
               "univariate")
})
