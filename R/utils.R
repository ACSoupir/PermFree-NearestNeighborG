# Internal helpers for normalising inputs ------------------------------------

#' Coerce user data to a spatstat point pattern
#'
#' @param data A \code{data.frame} with columns \code{x}, \code{y} (and
#'   optionally a mark column), or an object of class
#'   \code{\link[spatstat.geom]{ppp}}.
#' @param mark_col Name of the column holding marks when \code{data} is a
#'   data frame.
#' @return A \code{ppp} object.
#' @keywords internal
.prepare_ppp <- function(data, mark_col = NULL) {
  if (inherits(data, "ppp")) {
    return(spatstat.geom::as.ppp(data))
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame or a spatstat ppp object.", call. = FALSE)
  }
  if (!all(c("x", "y") %in% names(data))) {
    stop("'data' must contain columns 'x' and 'y'.", call. = FALSE)
  }
  if (nrow(data) < 3L) {
    stop("At least three points are required to build a convex hull window.", call. = FALSE)
  }
  win <- tryCatch(
    spatstat.geom::convexhull.xy(data$x, data$y),
    error = function(e) stop("Could not build a window from 'x'/'y': ",
                             conditionMessage(e), call. = FALSE)
  )
  mk <- NULL
  if (!is.null(mark_col)) {
    if (!mark_col %in% names(data)) {
      stop(sprintf("Mark column '%s' not found in 'data'.", mark_col), call. = FALSE)
    }
    mk <- as.factor(data[[mark_col]])
  }
  spatstat.geom::ppp(x = data$x, y = data$y, window = win, marks = mk)
}

#' Mark levels of a point pattern
#'
#' @param pp A \code{ppp} object.
#' @return Character vector of levels, or \code{NULL} when unmarked.
#' @keywords internal
.mark_levels <- function(pp) {
  m <- spatstat.geom::marks(pp)
  if (is.null(m)) return(NULL)
  levels(as.factor(m))
}

#' Subset a point pattern to one mark level
#'
#' @param pp A marked \code{ppp}.
#' @param lvl Mark level to keep.
#' @return A \code{ppp} containing only points with that mark.
#' @keywords internal
.subset_mark <- function(pp, lvl) {
  m <- spatstat.geom::marks(pp)
  if (is.null(m)) {
    stop("A mark level was requested but the point pattern is unmarked.", call. = FALSE)
  }
  keep <- as.factor(m) == lvl
  if (!any(keep)) {
    stop(sprintf("No points have mark '%s'. Available: %s", lvl,
                 paste(levels(as.factor(m)), collapse = ", ")), call. = FALSE)
  }
  pp[keep]
}

#' Resolve the radius grid used by spatstat
#'
#' @param radii User-supplied radii (or \code{NULL}).
#' @return The user radii, invisibly validated.
#' @keywords internal
.check_radii <- function(radii) {
  if (is.null(radii)) return(NULL)
  radii <- as.numeric(radii)
  if (length(radii) < 1L || anyNA(radii)) {
    stop("'radii' must be a non-empty numeric vector.", call. = FALSE)
  }
  if (radii[1] != 0) {
    stop("'radii' must start at 0 (a spatstat requirement).", call. = FALSE)
  }
  if (is.unsorted(radii, strictly = TRUE)) {
    stop("'radii' must be strictly increasing.", call. = FALSE)
  }
  radii
}

#' Pull the observed curve out of a spatstat fv object for one correction
#'
#' @param obs An \code{fv} object returned by Gest/Gcross.
#' @param correction One of "none", "rs", "han", "km".
#' @return Numeric vector of observed values.
#' @keywords internal
.obs_column <- function(obs, correction) {
  col <- switch(correction,
                none = "raw",
                rs   = "rs",
                han  = "han",
                km   = "km")
  if (!col %in% names(obs)) {
    stop(sprintf("spatstat did not return a '%s' column for correction = '%s'.",
                 col, correction), call. = FALSE)
  }
  as.numeric(obs[[col]])
}
