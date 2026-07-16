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

  # CHANGED (performance): O(1) name -> row-index lookup instead of scanning
  # rownames with which(rownames(boundary) == row_name) for every face,
  # which made the fill-in O(n^2) string comparisons.
  row_index <- stats::setNames(seq_along(filist), rownames(boundary))

  last_1 <- rep(NA, ncol(boundary))
  for(i in seq_along(filist)) {

    # create boundary matrix
    dt <- metadata[[i]]
    for (one_face in dt) {
      if (length(one_face) == 0) next  # vertices have an empty boundary
      row_idx <- row_index[paste(one_face, collapse = " ")]
      if (is.na(row_idx)) next         # face not in the filtration (as before)
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
#' @param last_1 The last 1 row index for each column in boundary matrix (after reduction).
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

  # CHANGED (bug fix): essential classes.
  #
  # Old condition:  !(i %in% pivot_owner)
  # The *values* of pivot_owner are killer COLUMN indices, so this only
  # excluded negative (killing) simplices. A positive simplex that was
  # already paired (its ROW is a pivot, i.e. !is.na(pivot_owner[i]))
  # slipped through and produced a spurious infinite bar. Example: for the
  # 4-point square, vertices 2..4 die at t = 1 but were also reported as
  # essential H0 classes.
  #
  # Correct definition: a simplex i is essential iff
  #   (a) its column reduces to zero  -> is.na(last_1[i])   (positive simplex)
  #   (b) its row is never a pivot    -> is.na(pivot_owner[i]) (never killed)
  for (i in seq_along(filist)) {
    if (is.na(last_1[i]) && is.na(pivot_owner[i]) &&
        !all(is.na(filist[[i]]$simplex))) {
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

#' Persistence pairs via sparse boundary reduction (for large filtrations)
#'
#' Same standard column-reduction algorithm as
#' \code{boundary_info()} + \code{extract_persistence_pairs()}, but each
#' column is stored as a sorted integer vector over GF(2) instead of a row of
#' a dense n x n matrix. Memory drops from O(n^2) to O(total number of
#' non-zeros), which is what makes filtrations with thousands to hundreds of
#' thousands of simplices (e.g. Flood or large Vietoris-Rips complexes)
#' feasible. Output is identical to the dense pipeline.
#'
#' @param filist Filtration list, each element includes simplex and time.
#' @return A data frame with columns: dim, birth, and death.
#'
#' @importFrom utils combn
#' @export
#' @examples
#' points <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1), ncol = 2)
#' filtration <- build_vr_filtration(points, eps_max=1.2)
#' pairs <- persistence_pairs(filtration)
persistence_pairs <- function(filist) {
  n <- length(filist)
  keys <- vapply(filist, function(x) paste(x$simplex, collapse = " "), "")
  index <- new.env(hash = TRUE, parent = emptyenv())
  for (i in seq_len(n)) assign(keys[i], i, envir = index)

  # sparse boundary columns: indices of the facets of each simplex
  cols <- vector("list", n)
  for (i in seq_len(n)) {
    s <- filist[[i]]$simplex
    if (length(s) == 1L) { cols[[i]] <- integer(0); next }
    fmat <- combn(s, length(s) - 1L)
    cols[[i]] <- sort(vapply(seq_len(ncol(fmat)), function(j)
      get(paste(fmat[, j], collapse = " "), envir = index), 0L))
  }

  # XOR of two sparse GF(2) columns = symmetric difference of index sets
  symdiff <- function(a, b) sort.int(c(a[!(a %in% b)], b[!(b %in% a)]))

  pivot_owner <- rep(NA_integer_, n)
  for (j in seq_len(n)) {
    col <- cols[[j]]
    repeat {
      if (length(col) == 0L) break
      piv <- col[length(col)]           # lowest 1 = last index
      owner <- pivot_owner[piv]
      if (is.na(owner)) { pivot_owner[piv] <- j; break }
      col <- symdiff(col, cols[[owner]])
    }
    cols[[j]] <- col
  }

  dims <- lengths(lapply(filist, `[[`, "simplex")) - 1L
  ts <- vapply(filist, `[[`, 0, "t")

  killed <- which(!is.na(pivot_owner))
  pairs <- data.frame(
    dim = dims[killed],
    birth = ts[killed],
    death = ts[pivot_owner[killed]]
  )
  # essential classes: zero column (positive simplex) that is never paired
  zero_cols <- which(vapply(cols, length, 0L) == 0L)
  essential <- zero_cols[is.na(pivot_owner[zero_cols])]
  if (length(essential)) {
    pairs <- rbind(pairs, data.frame(
      dim = dims[essential], birth = ts[essential], death = Inf))
  }
  rownames(pairs) <- NULL
  pairs[order(pairs$dim, pairs$birth), ]
}
