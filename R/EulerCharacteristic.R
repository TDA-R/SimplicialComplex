#' Compute the Euler characteristic \eqn{\chi} of a simplicial complex
#'
#' This function computes the Euler characteristic of a simplicial complex based on the alternating sum of Betti numbers.
#'
#' @param simplices A list of simplices (each a numeric vector).
#' @param eps Optional numerical tolerance to pass to \code{rankMatrix()}.
#'
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
#' @examples
#' simplices <- list(c(1,2,3), c(2,3,4), c(5))
#' euler_characteristic(simplices, eps = 1e-3)
#'
#' @seealso \code{\link{betti_number}}
#' @export
euler_characteristic <- function(
    simplices, eps
) {
  max_dim <- max(sapply(simplices, length))
  print(max_dim)
  #　χ = Σ (-1)^k * β_k
  euler_char <- sum(sapply(0:(max_dim), function(bound_dim) {
    (-1)^bound_dim * betti_number(simplices, bound_dim, eps)
  }))

  print(paste('Euler Characteristic', euler_char, sep = ': '))

  return(euler_char)
}
