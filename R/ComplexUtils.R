#' Build the clique (flag) complex of a graph on a fixed vertex set
#'
#' Shared by \code{\link{VietorisRipsComplex}}, \code{\link{CechComplex}} and \code{\link{WitnessComplex}}: each of them reduces to "connect vertices
#' that are close enough (in whatever sense that complex uses), then take the maximal cliques of that graph as the maximal simplices."
#'
#' @param n Number of vertices.
#' @param edges An integer matrix with two columns (or a length-\code{2k} vector, igraph style), one row per edge, 1-based vertex indices.
#' @return A list with \code{network} (an \code{igraph} object, the 1-skeleton) and \code{simplices} (a list of integer vectors, the
#'   vertex sets of the maximal cliques, each sorted).
#'
#' @keywords internal
build_clique_complex <- function(n, edges) {
  network <- igraph::make_empty_graph(n = n, directed = FALSE)
  if (length(edges) > 0) {
    edges_mat <- matrix(edges, ncol = 2)
    network <- igraph::add_edges(network, as.vector(t(edges_mat)))
  }
  cliques <- igraph::max_cliques(network)
  simplices <- lapply(cliques, function(clique) as.vector(sort(clique)))
  list(network = network, simplices = simplices)
}

#' Expand maximal simplices into a sorted filtration list
#'
#' Shared filtration-assembly step used by every \code{build_filtration()} method: take a set of maximal simplices, generate every face of every
#' dimension via \code{\link{faces}}, assign each face a filtration time via \code{scale_fn}, and sort by (time, dimension, lexicographic order).
#'
#' @param maximal_simplices A list of integer vectors (the maximal simplices).
#' @param scale_fn A function taking one simplex (integer vector) and returning its filtration time.
#' @param max_dimension Optional integer cap. If supplied, faces of dimension greater than \code{max_dimension} are never generated in the first
#'   place - this is a structural cap on \code{kmax}, evaluated BEFORE \code{\link{faces}} is called, not a post-hoc filter. Callers that need
#'   the persistence "+1 trick" (see \code{\link{restrict_filtration}}) are responsible for passing \code{max_dimension + 1} here, not
#'   \code{max_dimension} itself - this function does not know about that convention.
#' @return A filtration list: one \code{list(simplex, t)} per face, sorted by (t, dimension, lexicographic order).
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
#' @param query Optional second numeric matrix; if supplied, returns the \code{nrow(points)} x \code{nrow(query)} cross-distance matrix instead of
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
#' Used by \code{\link{CechComplex}}: the Cech complex includes a simplex \eqn{\sigma} at scale \eqn{\epsilon} exactly when the balls of radius
#' \eqn{\epsilon} centered at its vertices have a common point, which happens if and only if the minimum enclosing ball of \eqn{\sigma}'s vertices has
#' radius at most \eqn{\epsilon} (this is the standard reduction used e.g. by GUDHI's Cech complex; see Cavanna, Jahanseir and Sheehy (2017)).
#'
#' @param points A numeric matrix, one point per row (at least 1 row).
#' @return A list with \code{center} (numeric vector) and \code{radius}.
#'
#' @keywords internal
min_enclosing_ball <- function(points) {

  points <- as.matrix(points)
  d <- ncol(points)

  ball_of <- function(boundary) {
    # boundary: matrix of <= d+1 points that must lie on the returned ball.
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

  # shuffle for the algorithm's expected linear-time guarantee; use a fixed (non-random) rotation.
  n <- nrow(points)
  ord <- if (n <= 1) seq_len(n) else c(2:n, 1)
  welzl(points[ord, , drop = FALSE], points[0, , drop = FALSE])
}

#' Circumsphere of an affinely independent point set
#'
#' The unique sphere through \code{points} whose center lies in their affine hull, i.e. the minimal-radius sphere with all of \code{points} on its
#' boundary. Used by \code{\link{AlphaComplex}}/\code{\link{DelaunayComplex}} to compute alpha values (unlike \code{\link{min_enclosing_ball}}, which
#' minimizes radius over all enclosing balls.
#'
#' @param points A numeric matrix, one point per row (must be affinely independent.
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


#' Global GF(2) (Galois Field of order 2) boundary-matrix pivot reduction (Zomorodian-Carlsson)
#'
#' Shared reduction core used by \code{\link{persistence_pairs}}, \code{\link{flood_persistence}}, and \code{DiscreteMorse.R}'s internal \code{.dim1_triangle_pairing()}.
#' Given a filtration list (simplices in filtration order, each \code{list(simplex, t}), builds the sparse GF(2) boundary matrix,
#' column \code{j} holds the row indices (into \code{filist}) of the facets of \code{filist[[j]]},
#' and reduces it via standard pivot (low = highest surviving row index) elimination.
#'
#' @param filist A filtration list, each element must have \code{$simplex} (a vector of vertex ids).
#' @return A list with:
#'   \item{pivot_owner}{integer vector, length \code{length(filist)}. For row \code{i}, \code{pivot_owner[i]} is the column index \code{j}
#'     whose reduced column's pivot (lowest 1 / highest surviving row index) is \code{i}. \code{NA} if simplex \code{i} is never a pivot.}
#'   \item{cols}{list of the reduced sparse columns (integer row-index vectors), one per simplex. A zero-length reduced column marks a
#'     positive (creator) simplex; combined with \code{is.na(pivot_owner[i])} this identifies essential classes.}
#'
#' @keywords internal
.reduce_gf2_boundary <- function(filist) {

  n <- length(filist)
  keys <- vapply(filist, function(x) paste(x$simplex, collapse = " "), "")
  index <- new.env(hash = TRUE, parent = emptyenv())
  for (i in seq_len(n)) assign(keys[i], i, envir = index)

  # sparse boundary columns: indices of the facets of each simplex
  cols <- vector("list", n)
  for (i in seq_len(n)) {
    s <- filist[[i]]$simplex
    if (length(s) <= 1L) { cols[[i]] <- integer(0); next }
    fmat <- utils::combn(s, length(s) - 1L)
    cols[[i]] <- sort(vapply(seq_len(ncol(fmat)), function(j)
      get(paste(fmat[, j], collapse = " "), envir = index), 0L))
  }

  # XOR of two sparse GF(2) columns = symmetric difference of index sets
  symdiff <- function(a, b) sort.int(c(a[!(a %in% b)], b[!(b %in% a)]))

  pivot_owner <- rep(NA_integer_, n)
  for (j in seq_len(n)) {
    col <- cols[[j]]
    repeat {
      if (length(col) == 0L) break
      piv <- col[length(col)] # lowest 1 = last (highest) index
      owner <- pivot_owner[piv]
      if (is.na(owner)) { pivot_owner[piv] <- j; break }
      col <- symdiff(col, cols[[owner]])
    }
    cols[[j]] <- col
  }

  list(pivot_owner = pivot_owner, cols = cols)
}
