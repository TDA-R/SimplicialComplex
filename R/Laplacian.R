#' Persistent (combinatorial) Laplacian Delta_q^{X,Y}
#'
#' @param X_simplices,Y_simplices Lists of (maximal) simplices, same format used everywhere else in the package. X must be a subcomplex of Y.
#' @param q The dimension.
#' @return A list with the full Laplacian, its upper/down pieces, and the q-simplex basis (in row/column order) everything is expressed in.
persistent_laplacian <- function(X_simplices, Y_simplices, q) {

  simplex_key <- function(s) paste(sort(s), collapse = "-")
  simplex_label <- function(s) paste(sort(s), collapse = "")

  Xq <- faces(X_simplices, q)
  Yq <- faces(Y_simplices, q)
  Yq1 <- faces(Y_simplices, q + 1)

  nX <- length(Xq)
  if (nX == 0) stop(sprintf("X has no %d-simplices.", q))

  Xq_keys <- vapply(Xq, simplex_key, character(1)) # 0-1 1-2 2-3 0-3
  Yq_keys <- vapply(Yq, simplex_key, character(1)) # 0-1 0-2 1-2 0-3 2-3

  sel <- match(Xq_keys, Yq_keys) # 1 3 4 5 (no keys for 0-2)
  if (anyNA(sel)) {
    stop("X is not a subcomplex of Y at this dimension - some q-simplices of X don't appear in Y.")
  }
  extra <- setdiff(seq_along(Yq), sel) # only 2 (key for 0-2)

  # down: Depends only on X
  if (q == 0) {
    down <- matrix(0, nX, nX)
  } else {
    Dq_X <- boundary(X_simplices, q) # column order = faces(X, q) = Xq
    # [1,] -1  .  . -1
    # [2,]  1 -1  .  .
    # [3,]  .  1 -1  .
    # [4,]  .  .  1  1
    down <- as.matrix(Matrix::crossprod(Dq_X))
    #       [,1] [,2] [,3] [,4]
    # [1,]    2   -1    0    1
    # [2,]   -1    2   -1    0
    # [3,]    0   -1    2    1
    # [4,]    1    0    1    2
  }

  # upper: d_{q+1}^{X,Y} (d_{q+1}^{X,Y})*
  if (length(Yq1) == 0) {
    # Y has no (q+1)-simplices at all - C_{q+1}^{X,Y} is trivially {0}.
    upper <- matrix(0, nX, nX)
    dim_C <- 0
  } else {
    D <- as.matrix(boundary(Y_simplices, q + 1)) # rows = Yq, cols = Yq1
    #      [,1] [,2]
    # [1,]    1    0
    # [2,]   -1    1
    # [3,]    1    0
    # [4,]    0   -1
    # [5,]    0    1

    if (length(extra) == 0) {
      # Every q-simplex of Y is already in X, so nothing can leak out of X - all of C_{q+1}(Y) qualifies
      basis_C <- diag(length(Yq1))
    } else {
      # filter
      D_extra  <- D[extra, , drop = FALSE] # get the coef: -1 1
      basis_C  <- Null(t(D_extra)) # find the vertor ortho with this vec: 0.7071068 0.7071068
    }

    dim_C <- ncol(basis_C)
    if (dim_C == 0) {
      upper <- matrix(0, nX, nX)
    } else {
      # filter
      D_X <- D[sel, , drop = FALSE] # d_{q+1}^Y restricted to X's rows
      Mf  <- D_X %*% basis_C # d_{q+1}^{X,Y} in this basis

      # Because basis_C is orthonormal (P = I) and C_q(X) uses the canonical simplex basis (Q = I),
      # the adjoint is just a transpose no P^{-1} M^T Q correction needed.
      upper <- Mf %*% t(Mf)
    }
  }

  Delta <- upper + down
  labels <- vapply(Xq, simplex_label, character(1))
  dimnames(Delta) <- dimnames(upper) <- dimnames(down) <- list(labels, labels)

  list(laplacian = Delta, upper = upper, down = down, basis = Xq, labels = labels, dim_upper_chain = dim_C)
}


# Ordinary (non-persistent) Hodge Laplacian  L_k = B_k^T B_k + B_{k+1} B_{k+1}^T
#
# "Simplicial Attention Networks" This is the single-complex case:
# B_k and B_{k+1} both come from the SAME complex K, just two neighbouring dimensions of it - there is no X/Y split here.
#'
#' Ordinary Hodge Laplacian L_k(K) = B_k^T B_k + B_{k+1} B_{k+1}^T
#'
#' @param K_simplices A single complex (list of maximal simplices).
#' @param k The dimension.
#' @return A list with the full Laplacian, its down/up pieces, and the k-simplex basis (row/column order) everything is expressed in.
hodge_laplacian <- function(K_simplices, k) {

  simplex_label <- function(s) paste(sort(s), collapse = "")

  Kk <- faces(K_simplices, k)
  nK <- length(Kk)
  if (nK == 0) stop(sprintf("K has no %d-simplices.", k))

  # down: B_k^T B_k -> same complex K, dimension k boundary
  if (k == 0) {
    down <- matrix(0, nK, nK)
  } else {
    Bk   <- boundary(K_simplices, k) # column order = faces(K, k) = Kk
    down <- as.matrix(Matrix::crossprod(Bk))
  }

  # up: B_{k+1} B_{k+1}^T -> same complex K, dimension k+1 boundary
  Kk1 <- faces(K_simplices, k + 1)
  if (length(Kk1) == 0) {
    up <- matrix(0, nK, nK)
  } else {
    Bk1 <- as.matrix(boundary(K_simplices, k + 1)) # rows = Kk, cols = Kk1
    up  <- as.matrix(Matrix::tcrossprod(Bk1)) # Bk1 %*% t(Bk1)
  }

  Delta <- down + up
  labels <- vapply(Kk, simplex_label, character(1))
  dimnames(Delta) <- dimnames(down) <- dimnames(up) <- list(labels, labels)

  list(laplacian = Delta, up = up, down = down, basis = Kk, labels = labels)
}


# Basic graph Laplacian  L = D - A
#
# The classical graph-theory Laplacian - degree matrix minus adjacency
# matrix - built directly from a complex's 1-skeleton (vertices + edges).
#'
#' @param simplices A list of simplices; only the 0-simplices (vertices) and 1-simplices (edges) are used.
#' @return A list with the Laplacian L, the degree matrix D, the adjacency
#'   matrix A, and the vertex basis (row/column order).
graph_laplacian <- function(simplices) {

  simplex_label <- function(s) paste(sort(s), collapse = "")

  V <- faces(simplices, 0)
  E <- faces(simplices, 1)

  nV <- length(V)
  if (nV == 0) stop("Complex has no vertices.")

  labels <- vapply(V, simplex_label, character(1))
  idx <- setNames(seq_len(nV), labels)

  A <- matrix(0, nV, nV, dimnames = list(labels, labels))

  for (e in E) {
    endpoints <- sort(e)
    i <- idx[[simplex_label(endpoints[1])]]
    j <- idx[[simplex_label(endpoints[2])]]
    A[i, j] <- A[i, j] + 1
    A[j, i] <- A[j, i] + 1
  }

  D <- diag(rowSums(A), nV, nV)
  dimnames(D) <- list(labels, labels)

  L <- D - A

  list(laplacian = L, D = D, A = A, basis = V, labels = labels)
}
