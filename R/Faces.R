#' Generate all unique faces of a given dimension from simplices
#'
#' This function computes all unique \eqn{k}-faces (subsets of vertices) from the input simplices
#' and returns only those of the specified dimension.
#'
#' @param simplices A list of simplices (each a numeric vector).
#' @param target_dim The target dimension \eqn{k} for the faces (e.g., \code{0} for vertices, \code{1} for edges, etc.).
#'
#' @return A list of faces (each a numeric vector) of dimension \code{target_dim}.
#'
#' @details
#' The function generates all possible subsets (combinations) of each simplex, removes duplicates,
#' and filters them to only include those of length \code{target_dim + 1}.
#'
#' For example, a 2-simplex \code{c(1, 2, 3)} has three 1-dimensional faces (edges):
#' \code{c(1,2)}, \code{c(1,3)}, and \code{c(2,3)}, and three 0-dimensional faces (vertices):
#' \code{1}, \code{2}, and \code{3}.
#'
#' @importFrom gtools combinations
#' @export
faces <- function(
    simplices, target_dim
) {
  faceset <- list()
  facestr_set <- character()

  for (simplex in simplices) {
    node_counts <- length(simplex)

    for (r in node_counts:1) {
      combs <- combinations(n = node_counts, r = r, v = simplex)

      for (i in 1:nrow(combs)) {
        face <- sort(combs[i, ])
        face_str <- paste(face, collapse = "-")

        if (!(face_str %in% facestr_set)) {
          faceset[[length(faceset) + 1]] <- face
          facestr_set <- c(facestr_set, face_str)
        }
      }
    }
  }

  faces_in_dim <- Filter(function(face) length(face) == target_dim + 1, faceset)

  return(faces_in_dim)
}
