#' Construct a Cech Complex
#'
#' @param points A numeric matrix or data.frame with one point per row
#'   (columns are coordinates).
#' @param epsilon A non-negative numeric radius, used only to bound how many
#'   candidate simplices are generated (see Details) - not a promise that
#'   every returned simplex individually satisfies the Cech criterion at
#'   this scale.
#' @return A list with:
#' \describe{
#'   \item{network}{An \code{igraph} object: the 1-skeleton of the candidate
#'     graph (edges where distance \eqn{\le 2\epsilon}).}
#'   \item{simplices}{A list of integer vectors, each the vertex indices of
#'     a maximal clique of \code{network}.}
#' }
#'
#' @export
#' @examples
#' points <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1), ncol = 2)
#' cech_complex <- CechComplex(points, epsilon = 0.8)
CechComplex <- function(points, epsilon) {
  points <- as.matrix(points)
  n <- nrow(points)
  D <- pairwise_dist(points)

  edges <- which(upper.tri(D) & D <= 2 * epsilon, arr.ind = TRUE)
  build_clique_complex(n, edges)
}

#' @keywords internal
cech_scale_of_simplex <- function(points, simplex_idx) {
  if (length(simplex_idx) <= 1) return(0)
  min_enclosing_ball(points[simplex_idx, , drop = FALSE])$radius
}
