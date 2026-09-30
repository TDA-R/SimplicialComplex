#' Row-reduce (partial pivoting) a matrix
#'
#' @param M A numeric matrix.
#' @param tol Pivot values smaller than this (in absolute value) are treated as zero. Defaults to \code{1e-8}.
#' @return A list with \code{R} (the row-reduced matrix) and \code{pivots}
#'   (the column indices where a pivot was found, in row order).
#'
#' @details
#' Shared elimination step behind \code{\link{ker}} and \code{\link{im}}
#' (and, internally, the zigzag module's own linear solves): partial-pivoting
#' Gauss-Jordan elimination, stopping early once every row has a pivot.
#'
#' @export
gauss_jordan_eliminate <- function(M, tol = 1e-8) {
  M <- as.matrix(M)
  storage.mode(M) <- "double"
  nr <- nrow(M)
  nc <- ncol(M)
  pivots <- integer(0)
  r <- 1

  for (cc in seq_len(nc)) {
    if (r > nr) break
    col <- M[r:nr, cc]
    rel <- which.max(abs(col))
    piv_row <- rel + r - 1
    if (abs(M[piv_row, cc]) < tol) next
    if (piv_row != r) {
      tmp <- M[r, ]; M[r, ] <- M[piv_row, ]; M[piv_row, ] <- tmp
    }
    M[r, ] <- M[r, ] / M[r, cc]
    for (i in seq_len(nr)) {
      if (i != r && abs(M[i, cc]) > tol) M[i, ] <- M[i, ] - M[i, cc] * M[r, ]
    }
    pivots <- c(pivots, cc)
    r <- r + 1
  }
  list(R = M, pivots = pivots)
}

#' Compute a basis for the kernel (null space) of a matrix
#'
#' @param M A numeric matrix (or an object coercible to one, e.g. a sparse
#'   \code{Matrix}). Rows are the codomain, columns the domain: \eqn{M} is
#'   read as a linear map \eqn{M : \mathbb{R}^{\mathrm{ncol}(M)} \to
#'   \mathbb{R}^{\mathrm{nrow}(M)}}.
#' @param tol Pivoting tolerance passed to the Gauss-Jordan elimination used
#'   internally (see \code{\link{im}}, which shares the same elimination
#'   step). Defaults to \code{1e-8}.
#' @return A matrix with \code{ncol(M)} rows, one column per basis vector of
#'   \eqn{\ker(M)}. If the kernel is trivial (\eqn{\{0\}}), the result has
#'   \code{0} columns.
#'
#' @details
#' \eqn{\ker(M) = \{ x \in \mathbb{R}^{\mathrm{ncol}(M)} : Mx = 0 \}}. The
#' basis is obtained from Gauss-Jordan elimination of \eqn{M}: one basis
#' vector per free (non-pivot) column, in the usual parametric-solution
#' construction. If \code{M} has \code{0} rows (the zero map), every
#' standard basis vector of the domain is in the kernel, so the identity
#' matrix is returned.
#'
#' As with any basis, the specific vectors returned are not unique, they
#' depend on the pivoting order of the elimination, only the number of
#' columns (\eqn{\dim \ker(M)}) is an invariant of \code{M}.
#'
#' @export
ker <- function(M, tol = 1e-8) {
  M <- as.matrix(M)
  nr <- nrow(M)
  nc <- ncol(M)

  if (nc == 0) return(matrix(numeric(0), nrow = 0, ncol = 0))
  if (nr == 0) return(diag(nc))

  rr <- gauss_jordan_eliminate(M, tol = tol)
  piv <- rr$pivots
  R <- rr$R
  free <- setdiff(seq_len(nc), piv)

  if (length(free) == 0) return(matrix(numeric(0), nrow = nc, ncol = 0))

  basis <- vapply(free, function(fcol) {
    v <- numeric(nc)
    v[fcol] <- 1
    for (r in seq_along(piv)) v[piv[r]] <- -R[r, fcol]
    v
  }, numeric(nc))

  matrix(basis, nrow = nc)
}

#' Compute a basis for the image (column space) of a matrix
#'
#' @param M A numeric matrix (or an object coercible to one, e.g. a sparse
#'   \code{Matrix}).
#' @param tol Pivoting tolerance passed to the Gauss-Jordan elimination used
#'   internally. Defaults to \code{1e-8}.
#' @return A matrix with \code{nrow(M)} rows, one column per basis vector of
#'   \eqn{\mathrm{im}(M)}. If the image is trivial (\eqn{\{0\}}), the result
#'   has \code{0} columns.
#'
#' @details
#' \eqn{\mathrm{im}(M) = \{ Mx : x \in \mathbb{R}^{\mathrm{ncol}(M)} \}}, the
#' column space of \code{M}. The basis returned is a subset of \code{M}'s own
#' columns, specifically, the pivot columns found by Gauss-Jordan
#' elimination rather than synthetic linear combinations, so each basis
#' vector is directly interpretable as one of the original columns of
#' \code{M} (e.g. the boundary of one specific simplex).
#'
#' @export
im <- function(M, tol = 1e-8) {
  M <- as.matrix(M)
  nr <- nrow(M)
  nc <- ncol(M)

  if (nc == 0 || nr == 0) return(matrix(numeric(0), nrow = nr, ncol = 0))

  rr <- gauss_jordan_eliminate(M, tol = tol)
  piv <- rr$pivots

  if (length(piv) == 0) return(matrix(numeric(0), nrow = nr, ncol = 0))

  M[, piv, drop = FALSE]
}

#' Compute a homology basis from a cycle basis and a boundary basis
#'
#' @param Z A matrix whose columns form a basis of the cycle space
#'   \eqn{Z_k = \ker(\partial_k)} (typically the output of \code{\link{ker}}).
#' @param B A matrix whose columns form a basis of the boundary space
#'   \eqn{B_k = \mathrm{im}(\partial_{k+1})} (typically the output of
#'   \code{\link{im}}), living in the same ambient space as \code{Z} (i.e.
#'   \code{nrow(B) == nrow(Z)}). Defaults to an empty basis (\code{B_k = 0}).
#' @param tol Numerical tolerance passed to \code{\link{betti_number}}'s
#'   rank helper (\code{safe_rank}). Defaults to \code{NULL}, i.e.
#'   \code{Matrix::rankMatrix()}'s own default tolerance.
#' @return A matrix with \code{nrow(Z)} rows, one column per representative
#'   of a basis of the quotient \eqn{H_k = Z_k / B_k}. If \eqn{H_k = 0}, the
#'   result has \code{0} columns; \code{ncol()} of the result is the Betti
#'   number \eqn{\beta_k = \dim H_k}.
#'
#' @details
#' \eqn{H_k = \ker(\partial_k) / \mathrm{im}(\partial_{k+1})}: two cycles
#' represent the same homology class exactly when they differ by a boundary.
#' This walks the columns of \code{Z} in order, greedily keeping any column
#' that increases the rank of the span accumulated so far (starting from
#' \code{B}'s span) - i.e. any cycle that is not already a linear
#' combination of \code{B} and the cycles kept before it. The kept columns
#' are one representative chain per homology class.
#'
#' \code{\link{betti_number}} computes \eqn{\dim H_k} directly from ranks
#' (rank-nullity), without ever materializing \code{Z}, \code{B}, or a
#' homology basis - that is cheaper when only the count is needed. Use
#' \code{homology()} when the actual representative cycles matter (e.g. to
#' visualize or track a specific hole), not to recompute a Betti number.
#'
#' As with \code{\link{ker}} and \code{\link{im}}, the specific
#' representative chosen for each class depends on the order cycles in
#' \code{Z} are tested against the growing boundary span, and is not unique.
#'
#' @export
homology <- function(Z, B = matrix(numeric(0), nrow = nrow(Z), ncol = 0), tol = NULL) {
  Z <- as.matrix(Z)
  B <- as.matrix(B)

  if (ncol(Z) == 0) return(matrix(numeric(0), nrow = nrow(Z), ncol = 0))

  current <- B
  r_now <- if (ncol(current) == 0) 0 else as.numeric(safe_rank(current, tol))
  keep <- logical(ncol(Z))

  for (j in seq_len(ncol(Z))) {
    test <- cbind(current, Z[, j])
    r_new <- as.numeric(safe_rank(test, tol))
    if (r_new > r_now) {
      keep[j] <- TRUE
      current <- test
      r_now <- r_new
    }
  }

  Z[, keep, drop = FALSE]
}
