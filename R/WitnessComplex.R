#' Construct a Witness Complex
#'
#' @param points A numeric matrix or data.frame of witness points \eqn{S},
#'   one point per row.
#' @param landmarks Either an integer (the number of landmarks to select from
#'   \code{points} via \code{\link{generate_landmarks}}, i.e. Farthest-Point
#'   Sampling) or an integer vector of row indices into \code{points} to use
#'   directly as landmarks \eqn{L}.
#' @param epsilon A non-negative numeric scale.
#' @param nu Landmark-distance parameter; defaults to \code{1}.
#' @return A list with:
#' \describe{
#'   \item{network}{An \code{igraph} object: the 1-skeleton on the landmark
#'     vertex set.}
#'   \item{simplices}{A list of integer vectors (indices into \code{landmarks}
#'     / \code{landmark_indices}), each the vertices of a maximal simplex.}
#'   \item{landmarks}{The landmark coordinate matrix.}
#'   \item{landmark_indices}{Row indices into \code{points}, or \code{NULL}
#'     if \code{landmarks} was given as coordinates rather than indices.}
#'   \item{edge_birth}{A symmetric matrix giving, for every pair of
#'     landmarks, the exact scale at which their edge is witnessed (used by
#'     \code{\link{build_filtration}} to time every simplex exactly, rather
#'     than only checking membership at the single scale \code{epsilon}).}
#' }
#'
#' @export
#' @examples
#' points <- matrix(rnorm(200), ncol = 2)
#' witness_complex <- WitnessComplex(points, landmarks = 15, epsilon = 0.5)
WitnessComplex <- function(points, landmarks, epsilon, nu = 1) {
  points <- as.matrix(points)

  lm_idx <- NULL
  if (length(landmarks) == 1 && is.numeric(landmarks)) {
    lm <- generate_landmarks(points, as.integer(landmarks))
    L <- lm$landmarks
    lm_idx <- lm$indices
  } else if (is.numeric(landmarks) && all(landmarks == as.integer(landmarks))) {
    lm_idx <- as.integer(landmarks)
    L <- points[lm_idx, , drop = FALSE]
  } else {
    L <- as.matrix(landmarks)
  }

  m <- nrow(L)
  D <- pairwise_dist(points, L) # N witnesses x m landmarks

  m_nu <- if (nu <= 0) {
    rep(0, nrow(D))
  } else {
    apply(D, 1, function(row) sort(row)[min(nu, length(row))])
  }

  # exact birth scale of every landmark pair: the smallest epsilon at which
  # some witness satisfies max(d(x_i,s), d(x_j,s)) <= m_nu(s) + epsilon
  edge_birth <- matrix(Inf, m, m)
  if (m >= 2) {
    for (i in 1:(m - 1)) {
      for (j in (i + 1):m) {
        b <- min(pmax(D[, i], D[, j]) - m_nu)
        edge_birth[i, j] <- edge_birth[j, i] <- b
      }
    }
  }

  edges <- which(upper.tri(edge_birth) & edge_birth <= epsilon, arr.ind = TRUE)
  complex <- build_clique_complex(m, edges)

  list(network = complex$network, simplices = complex$simplices,
       landmarks = L, landmark_indices = lm_idx, edge_birth = edge_birth)
}

#' @keywords internal
witness_scale_of_simplex <- function(edge_birth, simplex_idx) {
  if (length(simplex_idx) <= 1) return(0)
  pairs <- utils::combn(simplex_idx, 2)
  max(edge_birth[t(pairs)])
}
