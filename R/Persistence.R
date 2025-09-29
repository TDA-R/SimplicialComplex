library(Matrix)

source('./R/Filtration.R')

points <- as.matrix(iris[1:5, 1:2])
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
  last_1 <- rep(NA, ncol(boundary))
  for(i in seq_along(F)) {

    # create boundary matrix
    dt <- metadata[[i]]
    for (one_face in dt) {
      row_name <- paste(one_face, collapse = " ")
      row_idx <- which(rownames(boundary) == row_name)
      boundary[row_idx, i] <- 1
    }

    # last boundary for reduction
    for_row <- boundary[, i]
    all_1 <- which(for_row == 1)
    if(length(all_1) > 0) {
      last_1[i] <- tail(all_1, 1)
    }
  }

  xor_vec <- function(a, b) {
    (a + b) %% 2
  }

  # record which row is used by which column
  pivot_owner <- rep(NA, nrow(boundary))
  for(j in seq_along(F)) {
    repeat {

      low_idx <- last_1[j] # last 1 row idx
      if(is.na(low_idx)) break
      owner_col <- pivot_owner[low_idx]

      if(is.na(owner_col)) {
        pivot_owner[low_idx] <- j
        break
      } else {
        # XOR owner_col to j
        boundary[, j] <- xor_vec(boundary[, j], boundary[, owner_col])
        # update last_1[j]
        for_row <- boundary[, j]
        all_1 <- which(for_row == 1)

        if(length(all_1) > 0) {
          last_1[j] <- tail(all_1, 1)
        } else {
          last_1[j] <- NA
        }
      }
    }
  }

  list(boundary = boundary, last_1 = last_1, pivot_owner = pivot_owner)

}

H0 <- persistent_pairs_for_dim(F, k = 0)

last_1 <- H0$last_1
boundary <- H0$boundary
pivot_owner <- H0$pivot_owner
