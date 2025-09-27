library(Matrix)

source('./R/Filtration.R')

persistent_pairs_for_dim <- function(
    F, k
    ) {
  idx_k <- which(sapply(F, function(x) length(x$simplex) == k + 1)) # k-simplices
  idx_km1 <- which(sapply(F, function(x) length(x$simplex) == k)) # (k-1)-simplices

  set_matrix <- matrix(numeric(0), 0, 2)
  set_name <- c("birth", "death")

  # if k=0, all 0-simplices are births with death = Inf
  if (k == 0) {
    if (length(idx_k) == 0) {
      return(list(pairs=set_matrix, inf=set_matrix))
    }
    births <- sapply(F[idx_k], function(x) x$t)
    inf <- cbind(birth=births, death=rep(Inf, length(births)))
    colnames(inf) <- set_name
    return(list(pairs=set_matrix, inf=inf))
  }

  # k>0
  if (length(idx_k) == 0 || length(idx_km1) == 0) {
    return(list(pairs=set_matrix, inf=set_matrix))
  }

  key <- function(s) paste(s, collapse = "-")
  km1_keys <- setNames(seq_along(idx_km1), sapply(F[idx_km1], function(x) key(x$simplex)))

  I <- integer(0)
  J <- integer(0)
  X <- integer(0)

  for (cj in seq_along(idx_k)) {
    s <- F[[idx_k[cj]]]$simplex
    for (i in seq_along(s)) {
      face <- s[-i] # (k-1)-face
      r <- km1_keys[ key(face) ]
      if (!is.na(r)) {
        I <- c(I, unname(r))
        J <- c(J, cj)
        X <- c(X, 1L)
      }
    }
  }

  B <- Matrix::sparseMatrix(i=I, j=J, x=X, dims=c(length(idx_km1), length(idx_k)))
  B2  <- B %% 2
  low <- rep(NA, ncol(B2))
  get_low <- function(j) {
    nz <- which(B2[, j] != 0)
    if (length(nz) == 0) NA else max(nz)
  }
  for (j in seq_len(ncol(B2))) {
    repeat {
      lj <- get_low(j)
      if (is.na(lj)) { low[j] <- NA; break }
      p  <- which(low == lj)
      if (length(p) == 0) { low[j] <- lj; break }
      B2[, j] <- (B2[, j] + B2[, p[1]]) %% 2
    }
  }

  # Pivot
  pairs_list <- list()
  for (j in seq_len(ncol(B2))) {
    r <- low[j]
    if (!is.na(r)) {
      pairs_list[[length(pairs_list) + 1]] <-
        c(F[[ idx_km1[r] ]]$t, F[[ idx_k[j] ]]$t)
    }
  }
  pairs <- if (length(pairs_list)) do.call(rbind, pairs_list) else set_matrix
  colnames(pairs) <- set_name

  # Inf
  used_rows  <- low[!is.na(low)]
  alive_rows <- setdiff(seq_len(length(idx_km1)), used_rows)
  inf_list <- lapply(alive_rows, function(r) c(F[[ idx_km1[r] ]]$t, Inf))
  inf <- if (length(inf_list)) do.call(rbind, inf_list) else set_matrix
  colnames(inf) <- set_name
  list(pairs=pairs, inf=inf)
}

points <- as.matrix(iris[1:10, 1:2])

F <- build_vr_filtration(points, eps_max=2)

H0 <- persistent_pairs_for_dim(F, k = 0)
H1 <- persistent_pairs_for_dim(F, k = 1)

H0
H1
