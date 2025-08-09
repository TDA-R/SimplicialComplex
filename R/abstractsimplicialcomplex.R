#' Compute Euler characteristic for an abstract simplicial complex
#'
#' @param simplices A list of simplices (each a numeric vector).
#' @param dimension Optional max dimension to compute up to.
#' @param tol Optional numerical tolerance to pass to \code{rankMatrix()}.
#'
#' @return The Euler characteristic \eqn{\chi}.
#'
#' @export
abstract_simplicial_complex <- function(
    simplices, dimension, tol=NULL
) {

  # check
  max_dim <- max(sapply(simplices, length)) -1
  # we need +-1 dimension to calculate boundary and betti number
  # so if the max simplices is 3 point, the max dimension is 2
  if (dimension < 0 || dimension > max_dim) {
    stop(sprintf("dimension must be between 0 and %d", max_dim))
  }

  # simplices
  simplices <- lapply(simplices, function(simplex) sort(simplex))

  # betti numbers
  euler_characteristic(simplices, tol)
}
