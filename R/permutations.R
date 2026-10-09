# Permutation reference ------------------------------------------------------

#' Column name holding the observed curve for a correction
#'
#' @param correction One of "none", "rs", "han", "km".
#' @return Character scalar.
#' @keywords internal
.obs_col_name <- function(correction) {
  switch(correction, none = "raw", rs = "rs", han = "han", km = "km")
}

#' Permutation null for the univariate G function
#'
#' Exact-by-simulation reference: repeatedly draw a random subset of size
#' \code{n} from all locations and recompute the observed statistic.
#'
#' @inheritParams exact_gest
#' @param nsim Number of permutations.
#' @return A list with \code{mean}, \code{var}, \code{r}, \code{nsim} and
#'   \code{correction}.
#' @seealso [exact_gest()]
#' @export
permute_gest <- function(data, radii = NULL, mark_col = NULL, marks_i = NULL,
                         correction = c("none", "rs", "han", "km"),
                         nsim = 999L, n_cores = 1L) {
  correction <- match.arg(correction)
  if (is.null(radii)) {
    stop("'radii' is required for permutation runs so that every replicate ",
         "uses the same radius grid.", call. = FALSE)
  }
  radii <- .check_radii(radii)
  pp <- .prepare_ppp(data, mark_col = mark_col)
  if (is.null(marks_i)) {
    pp_sub <- pp
  } else {
    pp_sub <- .subset_mark(pp, marks_i)
  }
  N <- spatstat.geom::npoints(pp)
  n_obs <- spatstat.geom::npoints(pp_sub)

  obs <- .gest_fv(pp_sub, radii = radii, correction = correction)
  rr <- as.numeric(obs$r)
  col <- .obs_col_name(correction)

  sims <- .lapply_cores(seq_len(as.integer(nsim)), function(b) {
    S <- sample.int(N, n_obs)
    g <- .gest_fv(pp[S], radii = radii, correction = correction)
    as.numeric(g[[col]])
  }, n_cores = n_cores)

  M <- do.call(rbind, sims)
  list(mean = colMeans(M), var = apply(M, 2L, stats::var),
       r = rr, nsim = as.integer(nsim), correction = correction)
}

#' Permutation null for bivariate cross G
#'
#' Exact-by-simulation reference: permute the mark labels among all locations,
#' keeping the type counts fixed, and recompute \eqn{G_{ij}}.
#'
#' @inheritParams exact_gcross
#' @param nsim Number of permutations.
#' @return A list with \code{mean}, \code{var}, \code{r}, \code{nsim} and
#'   \code{correction}.
#' @seealso [exact_gcross()]
#' @export
permute_gcross <- function(data, radii = NULL, mark_col = NULL,
                           marks_i, marks_j,
                           correction = c("none", "rs", "han", "km"),
                           nsim = 999L, n_cores = 1L) {
  correction <- match.arg(correction)
  if (is.null(radii)) {
    stop("'radii' is required for permutation runs so that every replicate ",
         "uses the same radius grid.", call. = FALSE)
  }
  radii <- .check_radii(radii)
  if (missing(marks_i) || missing(marks_j)) {
    stop("Both 'marks_i' and 'marks_j' are required.", call. = FALSE)
  }
  pp <- .prepare_ppp(data, mark_col = mark_col)

  if (identical(as.character(marks_i), as.character(marks_j))) {
    return(permute_gest(data = pp, radii = radii, marks_i = marks_i,
                        correction = correction, nsim = nsim,
                        n_cores = n_cores))
  }

  pp_i <- .subset_mark(pp, marks_i)
  n_i <- spatstat.geom::npoints(pp_i)
  pp_j <- .subset_mark(pp, marks_j)
  n_j <- spatstat.geom::npoints(pp_j)
  N <- spatstat.geom::npoints(pp)
  if (n_i + n_j > N) {
    stop("Disjoint label permutation requires n_i + n_j <= N.", call. = FALSE)
  }

  obs <- .gcross_fv(pp, marks_i, marks_j, radii = radii,
                    correction = correction)
  rr <- as.numeric(obs$r)
  col <- .obs_col_name(correction)
  W <- spatstat.geom::Window(pp)

  sims <- .lapply_cores(seq_len(as.integer(nsim)), function(b) {
    Si <- sample.int(N, n_i)
    rem <- setdiff(seq_len(N), Si)
    Sj <- sample(rem, n_j)
    mk <- rep("..other", N)
    mk[Si] <- "i"
    mk[Sj] <- "j"
    pp2 <- spatstat.geom::ppp(x = pp$x, y = pp$y, window = W,
                              marks = factor(mk))
    g <- spatstat.explore::Gcross(pp2, i = "i", j = "j", r = radii,
                                  correction = correction)
    as.numeric(g[[col]])
  }, n_cores = n_cores)

  M <- do.call(rbind, sims)
  list(mean = colMeans(M), var = apply(M, 2L, stats::var),
       r = rr, nsim = as.integer(nsim), correction = correction)
}
