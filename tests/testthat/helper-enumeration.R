# Shared test helpers --------------------------------------------------------

# Timing guard.  Every test must finish well under two minutes; individual
# tests are expected to stay far below this.
expect_under <- function(expr, secs = 60) {
  t0 <- Sys.time()
  force(expr)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (elapsed > secs) {
    stop(sprintf("Test exceeded its time budget: %.1f s (limit %.0f s).",
                 elapsed, secs), call. = FALSE)
  }
  invisible(elapsed)
}

# Exhaustive enumeration of the univariate null over all choose(N, n) subsets.
enum_uni <- function(D, radii, n) {
  N <- nrow(D)
  ci <- utils::combn(N, n)
  out <- matrix(NA_real_, ncol(ci), length(radii))
  for (c in seq_len(ncol(ci))) {
    S <- ci[, c]
    DS <- D[S, S, drop = FALSE]
    diag(DS) <- Inf
    out[c, ] <- vapply(radii, function(k) mean(rowSums(DS <= k) > 0), numeric(1))
  }
  out
}

# Exhaustive enumeration of the disjoint bivariate null.
enum_cross <- function(D, radii, ni, nj) {
  N <- nrow(D)
  ci <- utils::combn(N, ni)
  out <- NULL
  for (c in seq_len(ncol(ci))) {
    Si <- ci[, c]
    rem <- setdiff(seq_len(N), Si)
    cj <- utils::combn(length(rem), nj)
    for (c2 in seq_len(ncol(cj))) {
      Sj <- rem[cj[, c2]]
      out <- rbind(out, vapply(
        radii,
        function(k) mean(rowSums(D[Si, Sj, drop = FALSE] <= k) > 0),
        numeric(1)))
    }
  }
  out
}

pop_mean <- function(M) colMeans(M)

# Population variance over all enumerated configurations (denominator m, not
# the sample denominator used by stats::var).
pop_var <- function(M) {
  apply(M, 2L, function(z) mean((z - mean(z))^2))
}
