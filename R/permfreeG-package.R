#' permfreeG: Permutation-Free Exact Nearest-Neighbour G CSR
#'
#' Closed-form mean and variance of the sample-specific complete spatial
#' randomness (CSR) null distribution for the nearest-neighbour \eqn{G}
#' function, univariate (\code{\link[spatstat.explore]{Gest}}) and bivariate
#' cross (\code{\link[spatstat.explore]{Gcross}}).
#'
#' The null model permutes marks among a fixed set of spatial locations, so the
#' marked subset is a uniformly random subset of the observed point pattern.
#' This reproduces the permutation distribution exactly (up to numerical
#' precision) without enumerating permutations.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib permfreeG, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom utils combn
## usethis namespace: end
NULL
