#' Compute a CROCKER matrix for a time-varying point cloud
#'
#' CROCKER = "Contour Realization Of Computed k-dimensional hole Evolution in the Rips complex".
#' Treats the \eqn{k}-th Betti number as a function of two parameters at once.
#'
#' @param point_clouds A list of length \eqn{n}: \code{point_clouds[[i]]} is the point cloud observed at time step \eqn{i}.
#' @param dim The homology dimension \eqn{k} to track (\eqn{b_k}); e.g. \code{0} for connected components, \code{1} for loops.
#' @param method Complex type passed through to \code{\link{build_filtration}}.
#'   Defaults to \code{"VR"} (Vietoris-Rips), which is what the CROCKER plot was originally defined for, it only needs pairwise distances, not an
#'   ambient embedding, so it applies to any time-varying metric space. Other complex can be substituted,
#'   but the result is then a generalised Betti-surface rather than a "CROCKER plot" in the strict.
#' @param eps_max Maximum scale (epsilon) value. If \code{NULL} (default), it is set to the largest pairwise
#'   distance seen in any single frame.
#' @param n_eps Number of proximity values sampled uniformly from \code{0} to \code{eps_max} (the original CROCKER paper uses 50).
#' @param max_dimension Optional cap forwarded to \code{build_filtration()}/ \code{persistence_pairs()}. Defaults to \code{dim}.
#' @param n_cores Number of cores to use via \code{parallel::mclapply()} for the per-frame computation (see Details for why this parallelises
#'   cleanly). Defaults to \code{1} (sequential). Falls back to sequential on Windows, where \code{mclapply()} (fork-based) is unavailable.
#'
#' @return An object of class \code{"crocker"}: a list with
#' \describe{
#'   \item{matrix}{An \code{n_eps} x \code{n} numeric matrix; entry \code{[j, i]} is \eqn{b_k} of frame \eqn{i} at \code{eps_grid[j]}.}
#'   \item{eps_grid}{The sampled scale values (length \code{n_eps}).}
#'   \item{time}{Frame indices, \code{seq_along(point_clouds)}.}
#'   \item{dim}{The homology dimension tracked.}
#'   \item{long}{A tidy data frame with columns \code{t}, \code{epsilon}, \code{betti}, ready for \code{\link{plot_crocker}}.}
#' }
#'
#' @export
crocker <- function(
    point_clouds, dim, method = "VR", eps_max = NULL, n_eps = 50, max_dimension = NULL, n_cores = 1
    ) {

  if (!is.list(point_clouds) || length(point_clouds) == 0) {
    stop("point_clouds must be a non-empty list of point clouds, one per time step.")
  }

  point_clouds <- lapply(point_clouds, as.matrix)
  n_frames <- length(point_clouds)

  if (is.null(max_dimension)) max_dimension <- dim

  if (is.null(eps_max)) {
    eps_max <- max(vapply(point_clouds, function(pts) {
      if (nrow(pts) < 2) return(0)
      max(pairwise_dist(pts))
    }, numeric(1)))
    if (!is.finite(eps_max) || eps_max <= 0) {
      stop("Could not infer a positive eps_max from point_clouds; supply eps_max explicitly.")
    }
  }

  eps_grid <- seq(0, eps_max, length.out = n_eps)

  frame_column <- function(pts) {

    if (nrow(pts) == 0) return(rep(0, n_eps))
    fil <- build_filtration(pts, method = method, eps_max = eps_max,
                             max_dimension = max_dimension)
    pairs <- persistence_pairs(fil)
    pd <- pairs[pairs$dim == dim, , drop = FALSE]
    if (nrow(pd) == 0) return(rep(0, n_eps))
    vapply(eps_grid, function(e) sum(pd$birth <= e & pd$death > e), numeric(1))
  }

  use_parallel <- n_cores > 1 && .Platform$OS.type == "unix"
  columns <- if (use_parallel) {
    parallel::mclapply(point_clouds, frame_column, mc.cores = n_cores)
  } else {
    lapply(point_clouds, frame_column)
  }

  betti_mat <- do.call(cbind, columns)
  dimnames(betti_mat) <- NULL

  long <- data.frame(
    t = rep(seq_len(n_frames), each = n_eps),
    epsilon = rep(eps_grid, times = n_frames),
    betti = as.vector(betti_mat)
  )

  structure(
    list(matrix = betti_mat, eps_grid = eps_grid, time = seq_len(n_frames), dim = dim, long = long),
    class = "crocker"
  )
}

#' Plot a CROCKER matrix as a filled contour plot
#'
#' @param cr An object returned by \code{\link{crocker}}, or a data frame with columns \code{t}, \code{epsilon}, \code{betti} (e.g. \code{cr$long}).
#' @return A ggplot2 object: a tile plot of the Betti number over time and scale.
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_contour scale_fill_viridis_c labs theme_minimal theme .data
#' @export
plot_crocker <- function(cr) {

  df <- if (inherits(cr, "crocker")) cr$long else cr

  if (!all(c("t", "epsilon", "betti") %in% names(df))) {
    stop("plot_crocker() needs columns t, epsilon, betti (pass a crocker object or its $long data frame).")
  }

  fill_label <- if (inherits(cr, "crocker")) paste0("b", cr$dim) else "betti"

  contour_breaks <- seq(0.5, max(df$betti, 0) + 0.5, by = 1)

  ggplot(df, aes(x = .data$t, y = .data$epsilon)) +
    geom_tile(aes(fill = .data$betti)) +
    geom_contour(aes(z = .data$betti), color = "white", linewidth = 0.3, breaks = contour_breaks) +
    scale_fill_viridis_c(name = fill_label) +
    labs(title = "CROCKER plot", x = "Time", y = expression(paste("Proximity parameter ", epsilon))) +
    theme_minimal() +
    theme(legend.position = "right")
}
