#' Compute the Euler characteristic \eqn{\chi} of a simplicial complex
#'
#' @param simplices A list of simplices (each a numeric vector).
#' @param tol Optional numerical tolerance to pass to \code{rankMatrix()}.
#' @return An integer representing the Euler characteristic \eqn{\chi}.
#'
#' @details
#' The Euler characteristic is computed as:
#' \deqn{\chi = \sum_{k=0}^{k_{\max}} (-1)^k \beta_k}
#' where \eqn{\beta_k} is the \eqn{k}th Betti number, and \eqn{k_{\max}} is the highest dimension of any simplex in the complex.
#'
#' Interpretation of values:
#' \itemize{
#'   \item \eqn{\chi = 2}: Sphere-like surfaces
#'   \item \eqn{\chi = 1}: Disk-like spaces
#'   \item \eqn{\chi = 0}: Torus-like or circle-like spaces
#'   \item \eqn{\chi < 0}: Surfaces with multiple handles or genus
#' }
#'
#' @seealso \code{\link{betti_number}}
#' @export
#' @examples
#' simplices <- list(c(1, 2), c(3, 4), c(2, 1, 3), c(4, 2))
#' euler_characteristic(simplices, tol=0.1)
euler_characteristic <- function(
    simplices, tol
) {
  max_dim <- max(sapply(simplices, length))
  #　χ = Σ (-1)^k * β_k
  euler_char <- sum(sapply(0:(max_dim), function(bound_dim) {
    (-1)^bound_dim * betti_number(simplices, bound_dim, tol)
  }))

  return(euler_char)
}
