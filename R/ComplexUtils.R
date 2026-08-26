#' Build the clique (flag) complex of a graph on a fixed vertex set
#'
#' Shared by \code{\link{VietorisRipsComplex}}, \code{\link{CechComplex}} and
#' \code{\link{WitnessComplex}}: each of them reduces to "connect vertices
#' that are close enough (in whatever sense that complex uses), then take the
#' maximal cliques of that graph as the maximal simplices."
#'
#' @param n Number of vertices.
#' @param edges An integer matrix with two columns (or a length-\code{2k}
#'   vector, igraph style), one row per edge, 1-based vertex indices.
#' @return A list with \code{network} (an \code{igraph} object, the
#'   1-skeleton) and \code{simplices} (a list of integer vectors, the
#'   vertex sets of the maximal cliques, each sorted).
#'
#' @importFrom igraph graph.empty add_edges max_cliques
#' @keywords internal
build_clique_complex <- function(n, edges) {
  network <- graph.empty(n = n, directed = FALSE)
  if (length(edges) > 0) {
    edges_mat <- matrix(edges, ncol = 2)
    network <- add_edges(network, as.vector(t(edges_mat)))
  }
  cliques <- max_cliques(network)
  simplices <- lapply(cliques, function(clique) as.vector(sort(clique)))
  list(network = network, simplices = simplices)
}

#' Expand maximal simplices into a sorted filtration list
#'
#' Shared filtration-assembly step used by every \code{build_filtration()}
#' method: take a set of maximal simplices, generate every face of every
#' dimension via \code{\link{faces}}, assign each face a filtration time via
#' \code{scale_fn}, and sort by (time, dimension, lexicographic order) - the
#' same convention used throughout the package (see \code{\link{faces}} for
#' what "lexicographic order" means here).
#'
#' @param maximal_simplices A list of integer vectors (the maximal simplices).
#' @param scale_fn A function taking one simplex (integer vector) and
#'   returning its filtration time.
#' @param max_dimension Optional integer cap. If supplied, faces of dimension
#'   greater than \code{max_dimension} are never generated in the first
#'   place - this is a structural cap on \code{kmax}, evaluated BEFORE
#'   \code{\link{faces}} is called, not a post-hoc filter. Callers that need
#'   the persistence "+1 trick" (see \code{\link{restrict_filtration}}) are
#'   responsible for passing \code{max_dimension + 1} here, not
#'   \code{max_dimension} itself - this function does not know about that
#'   convention.
#' @return A filtration list: one \code{list(simplex = <integer vector>, t =
#'   <numeric>)} per face, sorted by (t, dimension, lexicographic order).
#'
#' @keywords internal
simplices_to_filtration <- function(maximal_simplices, scale_fn, max_dimension = NULL) {
  kmax <- max(sapply(maximal_simplices, length)) - 1
  if (!is.null(max_dimension)) {
    kmax <- min(kmax, max_dimension)
  }
  face <- lapply(0:kmax, function(k) faces(maximal_simplices, k))
  all_faces <- unlist(face, recursive = FALSE)

  simplex_info <- lapply(all_faces, function(s) list(simplex = s, t = scale_fn(s)))
  ord <- order(sapply(simplex_info, function(x) x[["t"]]),
               sapply(simplex_info, function(x) length(x$simplex)),
               sapply(simplex_info, function(x) paste(x$simplex, collapse = "-")))

  simplex_info[ord]
}

#' Pairwise Euclidean distance matrix
#'
#' @param points A numeric matrix, one point per row.
#' @param query Optional second numeric matrix; if supplied, returns the
#'   \code{nrow(points)} x \code{nrow(query)} cross-distance matrix instead of
#'   the full pairwise matrix of \code{points} with itself.
#' @return A numeric distance matrix.
#'
#' @keywords internal
pairwise_dist <- function(points, query = NULL) {
  a <- as.matrix(points)
  b <- if (is.null(query)) a else as.matrix(query)
  aa <- rowSums(a^2)
  bb <- rowSums(b^2)
  d2 <- outer(aa, bb, "+") - 2 * (a %*% t(b))
  sqrt(pmax(d2, 0))
}

#' Minimum enclosing ball of a finite point set (Welzl's algorithm)
#'
#' Used by \code{\link{CechComplex}}: the Cech complex includes a simplex
#' \eqn{\sigma} at scale \eqn{\epsilon} exactly when the balls of radius
#' \eqn{\epsilon} centered at its vertices have a common point, which happens
#' if and only if the minimum enclosing ball of \eqn{\sigma}'s vertices has
#' radius at most \eqn{\epsilon} (this is the standard reduction used e.g. by
#' GUDHI's Cech complex; see Cavanna, Jahanseir and Sheehy (2017)).
#'
#' Randomized, expected linear-time in the number of points; the point sets
#' passed in here are the vertices of a single simplex, so they are always
#' small.
#'
#' @param points A numeric matrix, one point per row (at least 1 row).
#' @return A list with \code{center} (numeric vector) and \code{radius}.
#'
#' @keywords internal
min_enclosing_ball <- function(points) {
  points <- as.matrix(points)
  d <- ncol(points)

  ball_of <- function(boundary) {
    # boundary: matrix of <= d+1 points that must lie ON the returned ball.
    m <- nrow(boundary)
    if (m == 0) return(list(center = rep(NA_real_, d), radius = -Inf))
    if (m == 1) return(list(center = boundary[1, ], radius = 0))
    b1 <- boundary[1, ]
    V <- sweep(boundary[-1, , drop = FALSE], 2, b1, "-")
    G <- V %*% t(V)
    q <- 0.5 * rowSums(V^2)
    lambda <- tryCatch(solve(G, q), error = function(e) {
      stop("Points are not in general position: cannot determine an enclosing ball.")
    })
    center <- b1 + as.vector(lambda %*% V)
    list(center = center, radius = sqrt(sum((center - b1)^2)))
  }

  in_ball <- function(p, ball, tol = 1e-9 * max(1, ball$radius)) {
    sqrt(sum((p - ball$center)^2)) <= ball$radius + tol
  }

  welzl <- function(pts, boundary) {
    if (nrow(pts) == 0 || nrow(boundary) == d + 1) {
      return(ball_of(boundary))
    }
    p <- pts[1, ]
    rest <- pts[-1, , drop = FALSE]
    ball <- welzl(rest, boundary)
    if (is.finite(ball$radius) && in_ball(p, ball)) return(ball)
    welzl(rest, rbind(boundary, p))
  }

  # shuffle for the algorithm's expected linear-time guarantee; use a fixed
  # (non-random) rotation instead of sample() so results are reproducible
  # without disturbing the caller's RNG state.
  n <- nrow(points)
  ord <- if (n <= 1) seq_len(n) else c(2:n, 1)
  welzl(points[ord, , drop = FALSE], points[0, , drop = FALSE])
}

#' Circumsphere of an affinely independent point set
#'
#' The unique sphere through \code{points} whose center lies in their affine
#' hull - i.e. the minimal-radius sphere with all of \code{points} on its
#' boundary. Used by \code{\link{AlphaComplex}}/\code{\link{DelaunayComplex}}
#' to compute alpha values (unlike \code{\link{min_enclosing_ball}}, which
#' minimizes radius over ALL enclosing balls, this fixes every point to lie
#' exactly on the sphere - the two coincide only for simplices whose
#' circumcenter is not "hidden" inside the simplex in a way that admits a
#' smaller covering ball).
#'
#' @param points A numeric matrix, one point per row (must be affinely
#'   independent, i.e. "general position" - see Otter et al. (2017)).
#' @return A list with \code{center} (numeric vector) and \code{radius}.
#'
#' @keywords internal
circumsphere <- function(points) {
  points <- as.matrix(points)
  m <- nrow(points)
  if (m == 1) return(list(center = points[1, ], radius = 0))

  b1 <- points[1, ]
  V <- sweep(points[-1, , drop = FALSE], 2, b1, "-")
  G <- V %*% t(V)
  q <- 0.5 * rowSums(V^2)
  lambda <- tryCatch(solve(G, q), error = function(e) {
    stop("Points are not in general position: cannot compute a circumsphere.")
  })
  center <- b1 + as.vector(lambda %*% V)
  list(center = center, radius = sqrt(sum((center - b1)^2)))
}
