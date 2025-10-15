#' Get the boundary matrix and its reduction information in matrix form
#'
#' @param filist Filtration list, each element includes simplex and time.
#' @return A list containing the boundary matrix, the last boundary row, and the pivot owner for persistence extraction.
#'
#' @importFrom utils combn tail
#' @export
#' @examples
#' points <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1), ncol = 2)
#' filtration <- build_vr_filtration(points, eps_max=1.2)
#' res <- boundary_info(filtration)
boundary_info <- function(filist) {

  boundary <- matrix(nrow = length(filist), ncol = length(filist), data = 0)
  name_vec = list()
  metadata = list()

  for (i in seq_along(filist)) {
    # add filist[[i]]$simplex to nam_vec
    name_vec[[i]] <- filist[[i]]$simplex
    # add info
    info_matrix <- combn(filist[[i]]$simplex, length(filist[[i]]$simplex) -1)
    metadata[[i]] <- lapply(seq_len(ncol(info_matrix)), function(j) info_matrix[, j])
  }

  rownames(boundary) <- colnames(boundary) <- sapply(name_vec, function(x) paste(x, collapse = " "))
  last_1 <- rep(NA, ncol(boundary))
  for(i in seq_along(filist)) {

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
  for(j in seq_along(filist)) {
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

  return(list(boundary = boundary, last_1 = last_1, pivot_owner = pivot_owner))

}

#' This function extracts the persistence from combining the boundary matrix and its filtration
#'
#' @param filist Filtration list, each element includes simplex and time.
#' @param last_1 The last 1 row index for each column in boundary matrix.
#' @param pivot_owner The column index owning the pivot row.
#' @return A data frame with columns: dimension, birth, and death.
#'
#' @export
#' @examples
#' points <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1), ncol = 2)
#' filtration <- build_vr_filtration(points, eps_max=1.2)
#' res <- boundary_info(filtration)
#' pairs <- extract_persistence_pairs(filtration, res$last_1, res$pivot_owner)
extract_persistence_pairs <- function(filist, last_1, pivot_owner) {
  pairs <- list()
  for (row in seq_along(pivot_owner)) {
    col <- pivot_owner[row]
    if (!is.na(col)) {
      # row: the simplex being killed, col: the simplex killing it
      birth_simplex <- filist[[row]]
      death_simplex <- filist[[col]]

      birth_time <- birth_simplex$t
      death_time <- death_simplex$t
      # the dimension of simplex that is killed
      dim <- length(birth_simplex$simplex) - 1

      pairs[[length(pairs)+1]] <- list(
        dim = dim,
        birth = birth_time,
        death = death_time,
        birth_simplex = birth_simplex$simplex,
        death_simplex = death_simplex$simplex
      )
    }
  }
  # not been killed (Inf)
  for (i in seq_along(filist)) {
    if (!(i %in% pivot_owner) && !all(is.na(filist[[i]]$simplex))) {
      birth_time <- filist[[i]]$t
      dim <- length(filist[[i]]$simplex) - 1
      pairs[[length(pairs)+1]] <- list(
        dim = dim,
        birth = birth_time,
        death = Inf,
        birth_simplex = filist[[i]]$simplex,
        death_simplex = NA
      )
    }
  }

  df <- data.frame(
    dim = sapply(pairs, function(x) x$dim),
    birth = sapply(pairs, function(x) x$birth),
    death = sapply(pairs, function(x) x$death)
  )
  return(df)
}



