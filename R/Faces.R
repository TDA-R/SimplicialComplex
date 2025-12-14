#' Generate all unique faces of a given dimension from simplices
#'
#' @param simplices A list of simplices (each a numeric vector).
#' @param target_dim The target dimension \eqn{k} for the faces (e.g., \code{0} for vertices, \code{1} for edges, etc.).
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
#' @examples
#' simplices <- list(c(1, 2), c(3, 4), c(2, 1, 3), c(4, 2))
#' faces(simplices, target_dim=0)
faces <- function(
    simplices, target_dim
) {

  faceset <- list()
  facestr_set <- character()

  target_len <- target_dim + 1

  for (simplex in simplices) {

    node_counts <- length(simplex)

    if (node_counts < target_len) {next}

    combs <- combinations(n = node_counts, r = target_len, v = simplex)

    for (i in 1:nrow(combs)) {
      face <- sort(combs[i, ])
      face_str <- paste(face, collapse = "-")

      if (!(face_str %in% facestr_set)) {
        faceset[[length(faceset) + 1]] <- face
        facestr_set <- c(facestr_set, face_str)
      }
    }
  }

  return(faceset)
}
