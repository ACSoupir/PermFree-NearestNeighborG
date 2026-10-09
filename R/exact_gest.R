# Univariate Gest ------------------------------------------------------------

#' Observed fv object for a correction
#'
#' @param X A \code{ppp}.
#' @param radii Radius vector (or \code{NULL} for the spatstat default grid).
#' @param correction One of "none", "rs", "han", "km".
#' @param rmax Optional maximum radius.
#' @return An \code{fv} object from [spatstat.explore::Gest()].
#' @keywords internal
.gest_fv <- function(X, radii = NULL, correction = "none", rmax = NULL) {
  args <- list(X = X, correction = correction)
  if (!is.null(radii)) args$r <- radii
  if (!is.null(rmax)) args$rmax <- rmax
  tryCatch(do.call(spatstat.explore::Gest, args),
           error = function(e) {
             stop("spatstat.explore::Gest() failed: ", conditionMessage(e),
                  call. = FALSE)
           })
}

#' Permutation-free exact CSR for the univariate nearest-neighbour G function
#'
#' Computes the observed \eqn{G(r)} for a marked subset together with the exact
#' mean and variance of that same statistic under random relabelling (a uniform
#' subset of the observed locations), removing any need to enumerate
#' permutations.
#'
#' @param data A \code{data.frame} with columns \code{x}, \code{y} (and
#'   optionally a mark column), or an object of class
#'   \code{\link[spatstat.geom]{ppp}}.
#' @param radii Numeric vector of radii, starting at 0.  \code{NULL} uses the
#'   spatstat default grid.
#' @param mark_col Name of the mark column when \code{data} is a data frame.
#' @param marks_i Mark level defining the subset.  \code{NULL} uses all points.
#' @param correction Edge correction: "none" (exact), "rs" (border, delta
#'   method from exact moments), or the mean-field "han"/"km" curves.
#' @param method \code{"exact"} (include the covariance terms) or
#'   \code{"approx"} (marginal variance only).
#' @param covariance Include off-diagonal covariance terms.  Ignored when
#'   \code{method = "approx"}.
#' @param n_cores Number of workers used for the C++ engine; currently only
#'   affects future parallel radius evaluation (reserved).
#' @param rmax Optional maximum radius passed to spatstat.
#' @return An object of class \code{nnG}: a data frame with columns
#'   \code{r}, \code{theo}, \code{obs}, \code{csr_mean}, \code{csr_var},
#'   \code{z} and \code{p_value}.
#' @seealso [exact_gcross()], [permute_gest()]
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(40), y = runif(40),
#'                 m = rep(c("A", "B"), each = 20))
#' res <- exact_gest(d, radii = seq(0, 0.5, by = 0.1), mark_col = "m",
#'                   marks_i = "A")
#' as.data.frame(res)
#' @export
exact_gest <- function(data, radii = NULL, mark_col = NULL, marks_i = NULL,
                       correction = c("none", "rs", "han", "km"),
                       method = c("exact", "approx"), covariance = TRUE,
                       n_cores = 1L, rmax = NULL) {
  correction <- match.arg(correction)
  method <- match.arg(method)
  radii <- .check_radii(radii)

  pp <- .prepare_ppp(data, mark_col = mark_col)
  if (is.null(marks_i)) {
    pp_sub <- pp
  } else {
    pp_sub <- .subset_mark(pp, marks_i)
  }
  n_obs <- spatstat.geom::npoints(pp_sub)
  if (n_obs < 2L) {
    stop("The marked subset must contain at least two points.", call. = FALSE)
  }

  obs <- .gest_fv(pp_sub, radii = radii, correction = correction, rmax = rmax)
  rr <- as.numeric(obs$r)

  cov_use <- isTRUE(covariance) && method == "exact"
  mom <- .analytic_moments(pp, rr, n_i = n_obs, correction = correction,
                           covariance = cov_use)

  df <- data.frame(r = rr,
                   theo = as.numeric(obs$theo),
                   obs = .obs_column(obs, correction),
                   csr_mean = as.numeric(mom$mean),
                   csr_var = as.numeric(mom$var))
  df$z <- ifelse(is.finite(df$csr_var) & df$csr_var > 0,
                 (df$obs - df$csr_mean) / sqrt(df$csr_var), NA_real_)
  df$p_value <- ifelse(is.na(df$z), NA_real_, 2 * stats::pnorm(-abs(df$z)))

  new_nnG(df, correction = correction,
          method = if (cov_use) "exact" else method,
          exact = correction == "none" && cov_use,
          N = spatstat.geom::npoints(pp), n_i = n_obs, statistic = "Gest")
}
