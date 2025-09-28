library(Matrix)

source('./R/Filtration.R')

points <- as.matrix(iris[1:10, 1:2])
F <- build_vr_filtration(points, eps_max=0.5)

persistent_pairs_for_dim <- function(
    F, k
    ) {

  boundary <- matrix(nrow = length(F), ncol = length(F), data = 0)
  name_vec = list()
  metadata = list()
  for (i in seq_along(F)) {
    # add F[[i]]$simplex to nam_vec
    name_vec[[i]] <- F[[i]]$simplex
    # add info
    info_matrix <- combn(F[[i]]$simplex, length(F[[i]]$simplex) -1)
    metadata[[i]] <- lapply(seq_len(ncol(info_matrix)), function(j) info_matrix[, j])
  }
  rownames(boundary) <- colnames(boundary) <- sapply(name_vec, function(x) paste(x, collapse = " "))

  for(i in seq_along(F)) {
    dt <- metadata[[i]]
    for (one_face in dt) {
      row_name <- paste(one_face, collapse = " ")
      row_idx <- which(rownames(boundary) == row_name)
      boundary[row_idx, i] <- 1
    }
  }

}

H0 <- persistent_pairs_for_dim(F, k = 0)
H0
