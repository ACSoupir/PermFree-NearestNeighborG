# Pure-R reference implementation -------------------------------------------
#
# Deliberately simple (explicit loops, no bitsets) so it can be read side by
# side with the derivation in vignette("derivation-univariate") and used to
# validate the fast C++ engine in the test suite.

#' Reference implementation of exact CSR moments (pure R)
#'
#' Computes the same quantities as [csr_uni_moments_cpp()] and
#' [csr_cross_moments_cpp()] using straightforward loops.  Intended for
#' validation and exposition, not production use.
#'
#' @param x,y Numeric vectors of point coordinates for the full pattern.
#' @param radii Radii at which to evaluate the moments.
#' @param n_i Subset size for univariate Gest, or type-i set size for bivariate
#'   Gcross.
#' @param n_j Type-j set size; \code{NULL} selects the univariate case.
#' @param covariance Include the off-diagonal covariance terms (exact variance)
#'   or only the marginal term.
#' @return A list with numeric vectors \code{mean} and \code{var}.
#' @keywords internal
csr_moments_r <- function(x, y, radii, n_i = NULL, n_j = NULL,
                          covariance = TRUE) {
  N <- length(x)
  cross <- !is.null(n_j)
  nsub_i <- as.integer(n_i)

  tg <- function(T, k) {
    out <- numeric(length(k))
    ok <- !is.na(k) & k >= 0 & k <= (length(T) - 1L)
    if (any(ok)) out[ok] <- T[k[ok] + 1L]
    out
  }

  D <- as.matrix(stats::dist(cbind(x, y)))
  T1 <- avoid_table(N - 1L, if (cross) as.integer(n_j) else nsub_i - 1L)
  T2 <- if (covariance) {
    avoid_table(N - 2L, if (cross) as.integer(n_j) else nsub_i - 2L)
  } else NULL
  c2 <- nsub_i * (nsub_i - 1) / (N * (N - 1))

  meanv <- varv <- numeric(length(radii))
  for (k in seq_along(radii)) {
    nb <- D <= radii[k]
    diag(nb) <- FALSE
    mp <- rowSums(nb)
    EX <- (nsub_i / N) * (1 - tg(T1, mp))
    meanv[k] <- mean(1 - tg(T1, mp))

    vv <- sum(EX * (1 - EX))
    if (covariance) {
      for (p in seq_len(N - 1L)) {
        for (q in (p + 1L):N) {
          if (!cross && nb[p, q]) {
            Eb <- c2
          } else if (cross) {
            mps <- mp[p] - as.numeric(nb[p, q])
            mqs <- mp[q] - as.numeric(nb[p, q])
            wc  <- sum(nb[p, ] & nb[q, ])
            Eb  <- c2 * (1 - tg(T2, mps) - tg(T2, mqs) +
                           tg(T2, mps + mqs - wc))
          } else {
            wc <- sum(nb[p, ] & nb[q, ])
            Eb <- c2 * (1 - tg(T2, mp[p]) - tg(T2, mp[q]) +
                          tg(T2, mp[p] + mp[q] - wc))
          }
          vv <- vv + 2 * (Eb - EX[p] * EX[q])
        }
      }
    }
    varv[k] <- max(vv / nsub_i^2, 0)
  }
  list(mean = meanv, var = varv)
}
