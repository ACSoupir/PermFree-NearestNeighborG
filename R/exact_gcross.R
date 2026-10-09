# Bivariate Gcross -----------------------------------------------------------

#' Observed fv object for a bivariate correction
#'
#' @param X A marked \code{ppp}.
#' @param i,j Mark levels.
#' @inheritParams .gest_fv
#' @return An \code{fv} object from [spatstat.explore::Gcross()].
#' @keywords internal
.gcross_fv <- function(X, i, j, radii = NULL, correction = "none",
                       rmax = NULL) {
  args <- list(X = X, i = i, j = j, correction = correction)
  if (!is.null(radii)) args$r <- radii
  if (!is.null(rmax)) args$rmax <- rmax
  tryCatch(do.call(spatstat.explore::Gcross, args),
           error = function(e) {
             stop("spatstat.explore::Gcross() failed: ", conditionMessage(e),
                  call. = FALSE)
           })
}

#' Permutation-free exact CSR for bivariate cross nearest-neighbour G
#'
#' Computes the observed \eqn{G_{ij}(r)} together with its exact mean and
#' variance under the disjoint label-permutation null: a random type-i set of
#' size \eqn{n_i} and, from the remaining locations, a random type-j set of
#' size \eqn{n_j}.
#'
#' @inheritParams exact_gest
#' @param marks_i,marks_j Mark levels defining the two types.  When they are
#'   equal the call is delegated to [exact_gest()] with \code{n = n_i}.
#' @return An object of class \code{nnG} (see [exact_gest()]).
#' @seealso [exact_gest()], [permute_gcross()]
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(60), y = runif(60),
#'                 m = sample(c("A", "B", "C"), 60, replace = TRUE))
#' res <- exact_gcross(d, radii = seq(0, 0.5, by = 0.1), mark_col = "m",
#'                     marks_i = "A", marks_j = "B")
#' as.data.frame(res)
#' @export
exact_gcross <- function(data, radii = NULL, mark_col = NULL,
                         marks_i, marks_j,
                         correction = c("none", "rs", "han", "km"),
                         method = c("exact", "approx"), covariance = TRUE,
                         n_cores = 1L, rmax = NULL) {
  correction <- match.arg(correction)
  method <- match.arg(method)
  radii <- .check_radii(radii)

  if (missing(marks_i) || missing(marks_j)) {
    stop("Both 'marks_i' and 'marks_j' are required.", call. = FALSE)
  }
  pp <- .prepare_ppp(data, mark_col = mark_col)

  if (identical(as.character(marks_i), as.character(marks_j))) {
    return(exact_gest(data = pp, radii = radii, marks_i = marks_i,
                      correction = correction, method = method,
                      covariance = covariance, n_cores = n_cores, rmax = rmax))
  }

  pp_i <- .subset_mark(pp, marks_i)
  pp_j <- .subset_mark(pp, marks_j)
  n_i <- spatstat.geom::npoints(pp_i)
  n_j <- spatstat.geom::npoints(pp_j)
  N <- spatstat.geom::npoints(pp)
  if (n_i < 2L) stop("Need at least two type-i points.", call. = FALSE)
  if (n_i + n_j > N) {
    stop("Disjoint label permutation requires n_i + n_j <= N.", call. = FALSE)
  }

  obs <- .gcross_fv(pp, marks_i, marks_j, radii = radii,
                    correction = correction, rmax = rmax)
  rr <- as.numeric(obs$r)

  cov_use <- isTRUE(covariance) && method == "exact"
  mom <- .analytic_moments(pp, rr, n_i = n_i, correction = correction,
                           covariance = cov_use, cross = TRUE, n_j = n_j)

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
          N = N, n_i = n_i, n_j = n_j, statistic = "Gcross")
}
