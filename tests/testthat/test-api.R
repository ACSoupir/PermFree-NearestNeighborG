test_that("data.frame and ppp inputs agree", {
  radii <- seq(0, 1.2, by = 0.3)
  d <- sim_nnG
  pp <- spatstat.geom::ppp(d$x, d$y,
                           window = spatstat.geom::convexhull.xy(d$x, d$y),
                           marks = d$m)
  a <- as.data.frame(exact_gest(d, radii = radii, mark_col = "m",
                                marks_i = "A"))
  b <- as.data.frame(exact_gest(pp, radii = radii, marks_i = "A"))
  expect_equal(a$csr_mean, b$csr_mean, tolerance = 1e-12)
  expect_equal(a$obs, b$obs, tolerance = 1e-12)
})

test_that("input validation is informative", {
  expect_error(exact_gest(sim_nnG, radii = seq(1, 2, by = 0.5),
                          mark_col = "m", marks_i = "A"),
               "start at 0")
  expect_error(exact_gest(sim_nnG, radii = seq(0, 1, by = 0.5),
                          mark_col = "m", marks_i = "Z"),
               "No points have mark")
  expect_error(exact_gest(sim_nnG, radii = seq(0, 1, by = 0.5),
                          mark_col = "nope", marks_i = "A"),
               "not found")
  expect_error(exact_gcross(sim_nnG, radii = seq(0, 1, by = 0.5),
                            mark_col = "m", marks_i = "A"),
               "marks_j")
  expect_error(permute_gest(sim_nnG, radii = NULL, mark_col = "m",
                            marks_i = "A", nsim = 10L),
               "radii")
})

test_that("nnG S3 methods work", {
  radii <- seq(0, 1.2, by = 0.3)
  res <- exact_gest(sim_nnG, radii = radii, mark_col = "m", marks_i = "A")
  expect_s3_class(res, "nnG")
  expect_output(print(res), "nnG")
  s <- summary(res)
  expect_s3_class(s, "summary.nnG")
  expect_output(print(s), "nnG summary")
  df <- as.data.frame(res)
  expect_true(all(c("r", "theo", "obs", "csr_mean", "csr_var", "z",
                    "p_value") %in% names(df)))
  grDevices::pdf(file = tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(res))
})

test_that("attributes record the analysis settings", {
  radii <- seq(0, 1.2, by = 0.3)
  res <- exact_gest(sim_nnG, radii = radii, mark_col = "m", marks_i = "A")
  expect_identical(res$correction, "none")
  expect_true(res$exact)
  expect_identical(res$N, 30L)
  expect_identical(res$n_i, 5L)

  app <- exact_gest(sim_nnG, radii = radii, mark_col = "m", marks_i = "A",
                    method = "approx")
  expect_false(app$exact)
})
