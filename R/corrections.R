# Edge-correction moments ----------------------------------------------------

#' Analytic CSR moments for one edge correction
#'
#' @param pp Full point pattern (all N locations).
#' @param radii Radii at which to evaluate the moments.
#' @param n_i Subset size (univariate) or type-i set size (bivariate).
#' @param correction One of "none", "rs", "han", "km".
#' @param covariance Include off-diagonal terms when exact.
#' @param cross Logical; bivariate Gcross?
#' @param n_j Type-j set size when \code{cross = TRUE}.
#' @return A list with numeric vectors \code{mean} and \code{var}
#'   (\code{var} is \code{NA} for the mean-field han/km curves).
#' @keywords internal
.analytic_moments <- function(pp, radii, n_i, correction = "none",
                              covariance = TRUE, cross = FALSE, n_j = NULL) {
  x <- pp$x
  y <- pp$y

  if (correction == "none") {
    return(if (cross) {
      csr_cross_moments_cpp(x, y, radii, n_i, as.integer(n_j), covariance)
    } else {
      csr_uni_moments_cpp(x, y, radii, n_i, covariance)
    })
  }

  if (correction == "rs") {
    bd <- spatstat.geom::bdist.points(pp)
    nj <- if (cross) as.integer(n_j) else -1L
    out <- rs_moments_cpp(x, y, bd, radii, n_i, nj)
    return(list(mean = out$mean, var = out$var))
  }

  if (correction %in% c("han", "km")) {
    if (cross) {
      stop("Analytic han/km moments are only implemented for univariate Gest; ",
           "use correction = 'none' or 'rs', or the permutation fallback.",
           call. = FALSE)
    }
    return(.meanfield_moments(pp, radii, n_i, correction))
  }

  stop("Unknown correction: ", correction, call. = FALSE)
}

#' Mean-field Hanisch / Kaplan-Meier CSR curve
#'
#' Uses the exact expectation of the uncensored event histogram and of the
#' observed-distance histogram, then assembles the corrected curve by plug-in.
#' The variance is not available analytically and is returned as \code{NA}.
#'
#' @inheritParams .analytic_moments
#' @return A list with \code{mean} and an all-\code{NA} \code{var}.
#' @keywords internal
.meanfield_moments <- function(pp, radii, n_subset, correction) {
  N <- spatstat.geom::npoints(pp)
  x <- pp$x
  y <- pp$y
  bd <- spatstat.geom::bdist.points(pp)
  Rr <- length(radii)

  T1 <- avoid_table(N - 1L, n_subset - 1L)
  p_hit <- function(k) {
    k2 <- as.integer(round(pmin(pmax(k, 0), N - 1L)))
    1 - T1[k2 + 1L]
  }

  # P(p in S, nnd <= min(r, b)) accumulated over radii
  mmin <- censored_neighbor_counts_cpp(x, y, bd, radii)
  cumP <- colSums((n_subset / N) * matrix(p_hit(mmin), nrow = N))
  Enco <- diff(c(0, cumP))                     # E[uncensored events per bin]

  # Distribution of o = min(nnd, b)
  rprev <- c(0, radii[-Rr])
  mprev <- censored_neighbor_counts_cpp(x, y, bd, rprev)
  Fm <- matrix(p_hit(mmin), nrow = N)
  Fp <- matrix(p_hit(mprev), nrow = N)
  mb <- counts_within_bdist_cpp(x, y, bd)
  Fb <- matrix(rep(p_hit(mb), Rr), nrow = N)
  inbin <- outer(bd, radii, function(bb, rr) as.numeric(bb <= rr)) -
           outer(bd, rprev, function(bb, rr) as.numeric(bb <= rr))
  Eobs <- colSums((n_subset / N) * (Fm - Fp + inbin * (1 - Fb)))

  if (correction == "km") {
    d <- rev(cumsum(rev(Eobs)))
    s <- ifelse(d > 0, 1 - Enco / d, 1)
    s[!is.finite(s)] <- 1
    curve <- 1 - cumprod(pmin(pmax(s, 0), 1))
  } else {
    W <- spatstat.geom::Window(pp)
    a <- spatstat.geom::eroded.areas(W, radii)
    good <- is.finite(a) & a > 0
    curve <- rep(NA_real_, Rr)
    if (any(good)) {
      incr <- numeric(Rr)
      incr[good] <- Enco[good] / a[good]
      last <- max(which(good))
      curve[seq_len(last)] <- cumsum(incr)[seq_len(last)]
      mx <- max(curve, na.rm = TRUE)
      if (is.finite(mx) && mx > 0) curve <- curve / mx
    }
  }
  list(mean = as.numeric(curve), var = rep(NA_real_, Rr))
}

#' Parallel lapply helper
#'
#' @param X List/vector to iterate over.
#' @param FUN Function applied to each element.
#' @param n_cores Number of workers (forked on unix, PSOCK otherwise).
#' @return A list.
#' @keywords internal
.lapply_cores <- function(X, FUN, n_cores = 1L) {
  if (n_cores > 1L && .Platform$OS.type == "unix") {
    parallel::mclapply(X, FUN, mc.cores = n_cores)
  } else if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::parLapply(cl, X, FUN)
  } else {
    lapply(X, FUN)
  }
}
