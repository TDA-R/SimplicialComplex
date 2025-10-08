#' Construct a Vietoris–Rips Complex (1-skeleton + maximal simplices)
#'
#' Builds the 1-skeleton graph by connecting pairs of points whose Euclidean
#' distance is strictly less than \eqn{\epsilon}, then extracts maximal cliques
#' as maximal simplices.
#'
#' @param points A numeric matrix or data.frame with one point per row (columns are coordinates).
#' @param epsilon A positive numeric threshold; connect points with distance < \eqn{\epsilon}.
#'
#' @return A list with:
#' \describe{
#'   \item{network}{An \code{igraph} object representing the 1-skeleton.}
#'   \item{simplices}{A list of integer vectors, each the vertex indices of a maximal simplex.}
#' }
#'
#' @details
#' The Vietoris–Rips complex at scale \eqn{\epsilon} includes a simplex for every
#' finite set of points with pairwise distances < \eqn{\epsilon}. This function
#' constructs the 1-skeleton (edges only) and then uses maximal cliques in that
#' graph as the maximal simplices.
#'
#' @importFrom igraph graph.empty add_edges max_cliques
#' @export
#' @examples
#' points <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1), ncol = 2)
#' epsilon <- 1.5
#' vr_complex <- VietorisRipsComplex(points, epsilon)
VietorisRipsComplex <- function(
    points, epsilon
    ) {

  labels <- 1:nrow(points)

  # Init graph (k-skeleton)
  network <- graph.empty(n = length(labels), directed = FALSE)

  for (i in 1:(length(labels) - 1)) {
    for (j in (i + 1):length(labels)) {
      dist <- sqrt(sum((points[i, ] - points[j, ])^2))
      if (dist <= epsilon) {
        # Add edge
        network <- add_edges(network, c(i, j))
      }
    }
  }
  cliques <- max_cliques(network)
  print(cliques)
  simplices <- lapply(cliques, function(clique) {
    simplex <- sort(clique)
    return(as.vector(simplex))
  })

  return(list(network = network, simplices = simplices))
}
