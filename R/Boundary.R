#' Compute the boundary operator for a simplicial complex
#'
#' This function returns the sparse matrix representation of the boundary operator \eqn{\partial_{k}}.
#'
#' @param simplices A list of simplices (each a numeric vector).
#' @param bound_dim The dimension k of the boundary operator \eqn{\partial_{k}}.
#'
#' @return A sparse matrix representing \eqn{\partial_{k}}.
#'
#' @details
#' \deqn{\partial_k \sigma = \sum_i (-1)^i [v_0 v_1 \ldots \hat{v}_i \ldots v_k]}
#'
#' @importFrom Matrix sparseMatrix Matrix
#' @importFrom stats setNames
#' @export
#' @examples
#' simplices <- list(c(1, 2), c(3, 4), c(2, 1, 3), c(4, 2))
#' boundary(simplices, 0)
boundary <- function(
    simplices, bound_dim
) {
  simplices_dim <- faces(simplices, target_dim=bound_dim)
  simplices_lower_dim <- faces(simplices, target_dim=(bound_dim-1))


  # No simplices at all in this dimension, return 0x0 matrix
  if (length(simplices_dim) == 0) {
    message(sprintf("No %d-simplices found, returning 0x0 matrix.", bound_dim))
    return(sparseMatrix(i = integer(0), j = integer(0), dims = c(0, 0)))
  }

  # No lower-dimensional simplices, return identity matrix
  if (length(simplices_lower_dim) == 0) {
    n <- length(simplices_dim)
    message(sprintf("No %d-dimensional lower simplices, returning identity matrix.", bound_dim - 1))
    i <- rep(1, n) # row
    j <- 1:n # column
    x <- rep(1, n) # value
    Sparse <- sparseMatrix(i=i, j=j, x=x)

  }
  else {

    # add index to simplices
    keys_dim <- sapply(simplices_dim, function(x) paste(sort(x), collapse= "-"))
    dct_dim <- setNames(seq_along(simplices_dim), keys_dim)

    keys_lower_dim <- sapply(simplices_lower_dim, function(x) paste(sort(x), collapse= "-"))
    dct_lower_dim <- setNames(seq_along(simplices_lower_dim), keys_lower_dim)

    Sparse <- Matrix(
      0,
      nrow = length(simplices_lower_dim),
      ncol = length(simplices_dim),
      sparse = TRUE
    )

    for (simplex_dim in simplices_dim) {
      j_key <- paste(simplex_dim, collapse = "-")
      j <- dct_dim[[j_key]]

      i_entries <- list()

      for (face in seq_along(simplex_dim)) {

        simplex_lower <- simplex_dim[-face]
        i_key <- paste(simplex_lower, collapse = "-")
        i <- dct_lower_dim[[i_key]]
        i_entries[[length(i_entries) + 1]] <- list(i = i, key = i_key)

        # Fill boundary matrix
        if (!is.na(i) && !is.na(j)) {
          Sparse[i, j] <- ifelse(face %% 2 == 0, -1, 1)
        }
      }

      # Print i and j only once per simplex
      i_message <- paste(sapply(i_entries, function(e) {
        sprintf("value: %d, key: %s", e$i, e$key)
      }), collapse = " I ")
      message("Lower dimension Simplices - ", i_message)
      message(sprintf("Simplices - value: %d, key: %s", j, j_key))
    }

  }
  return(Sparse)
}
