# Binomial / avoidance-table helpers ---------------------------------------
#
# The exact moments only ever need the ratio
#
#   T_{B,s}(k) = choose(B - k, s) / choose(B, s),   k = 0, ..., B
#
# which is the probability that a uniformly chosen s-subset of a set of size B
# avoids k designated elements.  Evaluating it from a precomputed table (once per
# analysis) instead of one lchoose() call per point pair is the single biggest
# speed-up over the original implementation.

#' Log binomial coefficient in log space
#'
#' \code{log_comb(n, k) = lgamma(n + 1) - lgamma(k + 1) -
#' lgamma(n - k + 1)}. Used to avoid overflow for large \code{n}.
#'
#' @param n,k Numeric vectors (recycled).
#' @return Numeric vector of log-binomial coefficients.
#' @keywords internal
log_comb <- function(n, k) {
  lgamma(n + 1) - (lgamma(k + 1) + lgamma(n - k + 1))
}

#' Avoidance table
#'
#' \code{avoid_table(B, s)[k + 1] = choose(B - k, s) / choose(B, s)}, with
#' value zero whenever \code{B - k < s}.
#'
#' @param B,s Non-negative integers.
#' @return Numeric vector of length \code{B + 1}.
#' @keywords internal
avoid_table <- function(B, s) {
  B <- as.integer(B)
  s <- as.integer(s)
  stopifnot(length(B) == 1L, length(s) == 1L, B >= 0L)
  k <- 0:B
  out <- numeric(length(k))
  ok <- (B - k) >= s & s >= 0L
  if (any(ok)) {
    out[ok] <- exp(lchoose(B - k[ok], s) - lchoose(B, s))
  }
  out
}
