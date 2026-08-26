#' Compute the persistence landscape of a persistence diagram
#'
#' Converts a persistence diagram data frame (the direct output of
#' \code{extract_persistence_pairs()}, \code{persistence_pairs()}, or
#' \code{flood_persistence()}) into its persistence landscape: a collection of
#' continuous, piecewise-linear functions \eqn{\lambda(k, \cdot)}, obtained by
#' overlaying the "tent" function of every birth-death pair and, at each time
#' \code{t}, taking the \eqn{k}th largest tent value (Bubenik, 2015; Chazal
#' and Michel, 2021, Section 5.4).
#'
#' @param df A persistence diagram data frame with columns \code{dim},
#'   \code{birth}, \code{death} (e.g. the output of
#'   \code{extract_persistence_pairs()}).
#' @param dimension Homological dimension to extract the landscape for.
#' @param k_max Number of landscape levels \eqn{\lambda(1,\cdot), \dots,
#'   \lambda(k_{max},\cdot)} to return. Defaults to all levels supported by
#'   the diagram (the number of finite birth-death pairs in \code{dimension});
#'   levels beyond that are identically zero.
#' @param resolution Number of grid points used to discretize each landscape
#'   function.
#' @param t_range Optional \code{c(t_min, t_max)} grid range. Defaults to
#'   \code{c(min(birth), max(death))} over the selected pairs.
#' @return A data frame with columns \code{t}, \code{k}, \code{value}: the
#'   discretized landscape functions, one row per (grid point, level) pair.
#'
#' @export
#' @examples
#' points <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1), ncol = 2)
#' filtration <- build_filtration(points, method = "VR", eps_max = 1.2)
#' res <- boundary_info(filtration)
#' pairs <- extract_persistence_pairs(filtration, res$last_1, res$pivot_owner)
#' landscape <- persistence_landscape(pairs, dimension = 0)
persistence_landscape <- function(
    df, dimension = 0, k_max = NULL, resolution = 500, t_range = NULL
) {

  sub <- df[df$dim == dimension, , drop = FALSE]

  n_inf <- sum(is.infinite(sub$death))
  if (n_inf > 0) {
    message(sprintf(
      "Dropping %d essential (death = Inf) pair(s) in dimension %d; persistence landscapes are only defined for finite bars.",
      n_inf, dimension
    ))
    sub <- sub[is.finite(sub$death), , drop = FALSE]
  }

  if (nrow(sub) == 0) {
    stop(sprintf("No finite birth-death pairs found in dimension %d.", dimension))
  }

  b <- sub$birth
  d <- sub$death
  n <- length(b)

  if (is.null(t_range)) {
    t_min <- min(b)
    t_max <- max(d)
  } else {
    t_min <- t_range[1]
    t_max <- t_range[2]
  }
  t_seq <- seq(t_min, t_max, length.out = resolution)

  if (is.null(k_max)) k_max <- n
  k_max <- min(k_max, n)

  # tent function values: tent[i, j] = Lambda_{(b_i, d_i)}(t_seq[j])
  tent <- matrix(0, nrow = n, ncol = resolution)
  for (i in seq_len(n)) {
    tent[i, ] <- pmax(0, pmin(t_seq - b[i], d[i] - t_seq))
  }

  # landscape[k, j] = kth largest tent value at t_seq[j]
  landscape <- matrix(0, nrow = k_max, ncol = resolution)
  for (j in seq_len(resolution)) {
    landscape[, j] <- sort(tent[, j], decreasing = TRUE)[seq_len(k_max)]
  }

  data.frame(
    t = rep(t_seq, each = k_max),
    k = rep(seq_len(k_max), times = resolution),
    value = as.vector(landscape)
  )
}
