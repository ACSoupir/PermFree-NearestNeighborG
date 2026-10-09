# nnG result class -----------------------------------------------------------

#' Construct an \code{nnG} object
#'
#' @param df Data frame with columns \code{r}, \code{theo}, \code{obs},
#'   \code{csr_mean}, \code{csr_var}, \code{z}, \code{p_value}.
#' @param correction,method Settings used.
#' @param exact Logical; was the analytic result fully exact?
#' @param N,n_i,n_j Point counts.
#' @param statistic \code{"Gest"} or \code{"Gcross"}.
#' @return An object of class \code{nnG}.
#' @keywords internal
new_nnG <- function(df, correction, method, exact, N,
                    n_i = NA_integer_, n_j = NA_integer_,
                    statistic = "Gest") {
  structure(list(data = df,
                 correction = correction,
                 method = method,
                 exact = isTRUE(exact),
                 N = as.integer(N),
                 n_i = as.integer(n_i),
                 n_j = if (is.null(n_j)) NA_integer_ else as.integer(n_j),
                 statistic = statistic),
            class = "nnG")
}

#' Methods for nnG objects
#'
#' @param x,object An \code{nnG} object.
#' @param y Ignored; present for compatibility with the plot generic.
#' @param row.names,optional,... Further arguments (ignored or passed on).
#' @return \code{as.data.frame} returns the underlying data frame,
#'   \code{summary} a list of class \code{summary.nnG}, and
#'   \code{print}, \code{plot} return their input invisibly.
#' @name nnG-methods
NULL

#' @rdname nnG-methods
#' @export
as.data.frame.nnG <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$data
}

#' @rdname nnG-methods
#' @export
print.nnG <- function(x, ...) {
  cat(sprintf("<nnG> %s | correction = '%s' | method = '%s'%s\n",
              x$statistic, x$correction, x$method,
              if (x$exact) " | exact" else ""))
  cat(sprintf("N = %d, n_i = %s%s\n", x$N,
              if (is.na(x$n_i)) "-" else as.character(x$n_i),
              if (!is.na(x$n_j)) paste0(", n_j = ", x$n_j) else ""))
  df <- x$data
  show <- seq_len(min(6L, nrow(df)))
  print(df[show, , drop = FALSE], row.names = FALSE)
  if (nrow(df) > 6L) cat("...", nrow(df), "rows total\n")
  invisible(x)
}

#' @rdname nnG-methods
#' @export
summary.nnG <- function(object, ...) {
  df <- object$data
  out <- list(statistic = object$statistic,
              correction = object$correction,
              method = object$method,
              exact = object$exact,
              N = object$N, n_i = object$n_i, n_j = object$n_j,
              mean_obs = mean(df$obs, na.rm = TRUE),
              max_z = if (all(is.na(df$z))) NA_real_ else
                df$z[which.max(abs(df$z))],
              min_p = if (all(is.na(df$p_value))) NA_real_ else
                min(df$p_value, na.rm = TRUE))
  class(out) <- "summary.nnG"
  out
}

#' @export
print.summary.nnG <- function(x, ...) {
  cat(sprintf("<nnG summary> %s | correction = '%s'%s\n", x$statistic,
              x$correction, if (x$exact) " | exact" else ""))
  cat(sprintf("N = %d, n_i = %s%s\n", x$N,
              if (is.na(x$n_i)) "-" else as.character(x$n_i),
              if (!is.na(x$n_j)) paste0(", n_j = ", x$n_j) else ""))
  cat(sprintf("mean observed G(r) = %.4f\n", x$mean_obs))
  if (!is.na(x$max_z)) cat(sprintf("largest |z| = %.3f (p = %g)\n", abs(x$max_z), x$min_p))
  invisible(x)
}

#' @rdname nnG-methods
#' @export
plot.nnG <- function(x, y = NULL, ...) {
  df <- x$data
  ylim <- range(c(df$obs, df$csr_mean, df$theo), na.rm = TRUE)
  graphics::plot(df$r, df$obs, type = "l", ylim = ylim, xlab = "r",
                 ylab = if (identical(x$statistic, "Gest")) "G(r)" else "Gcross(r)",
                 main = sprintf("%s CSR (correction = '%s')", x$statistic,
                                x$correction), ...)
  graphics::lines(df$r, df$theo, col = "grey40", lty = 2)
  graphics::lines(df$r, df$csr_mean, col = "red")
  if (any(is.finite(df$csr_var))) {
    sd <- sqrt(pmax(df$csr_var, 0))
    graphics::lines(df$r, df$csr_mean - sd, col = "red", lty = 3)
    graphics::lines(df$r, df$csr_mean + sd, col = "red", lty = 3)
  }
  graphics::legend("bottomright",
                   legend = c("observed", "theoretical CSR", "exact CSR mean"),
                   col = c("black", "grey40", "red"), lty = c(1, 2, 1), bty = "n")
  invisible(x)
}