#' Safely compute the rank of a sparse matrix
#'
#' This helper function wraps \code{Matrix::rankMatrix()} to safely handle empty matrices (i.e., with 0 rows or columns).
#'
#' @param simplices A list of simplices representing the simplicial complex.
#' @param bound_dim The dimension of the boundary to compute the Betti number for.
#' @param tol Optional numerical tolerance to pass to \code{rankMatrix()}.
#'
#' @return An integer representing the rank of the matrix.
#'
#' @keywords internal
#' @importFrom Matrix rankMatrix
betti_number <- function(
    simplices, bound_dim, tol = NULL
) {
  # Compute boundary matrices: \partial_k and \partial_{k+1}
  partial_i  <- boundary(simplices, bound_dim)
  partial_i1 <- boundary(simplices, bound_dim + 1)

  # Compute rank of \partial_{k}
  partial_i_rank <- if (bound_dim == 0) {
    0
  } else {
    safe_rank(partial_i, tol)
  }

  # Compute rank of \partial_{k+1}
  partial_i1_rank <- safe_rank(partial_i1, tol)

  # Compute Betti number: \beta_k = dim(C_k) - rank(\partial_k) - rank(\partial_{k+1})
  betti_num <- as.numeric(ncol(partial_i) - partial_i_rank - partial_i1_rank)

  print(paste("Betti", bound_dim, "=", betti_num, sep = ""))
  return(betti_num)
}

safe_rank <- function(mat, tol = NULL) {
  if (prod(dim(mat)) == 0) {
    return(0)
  }
  if (is.null(tol)) {
    return(rankMatrix(mat))
  } else {
    return(rankMatrix(mat, tol = tol))
  }
}

