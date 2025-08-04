#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

library(gtools)

abstract_simplicial_complex <- function(
  simplices, dimension
  ) {
  # simplices
  simplices <- lapply(simplices, function(simplex) sort(simplex))
  print(simplices)

  # faces
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

  faces_dim <- Filter(function(face) length(face) == dimension + 1, faceset)

  print("faces")
  print(facestr_set)
  print(faces_dim)

  # boundary

  # betti numbers

  # euler characteristic
}

simplices <- list(c(2, 1, 3), c(4, 2), c(5))
abstract_simplicial_complex(simplices, 1)
