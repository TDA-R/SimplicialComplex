#' @keywords internal
#' Example:
#'   .zz_simplex_key(c(2, 1, 3))
#'   #> [1] "1-2-3"
.zz_simplex_key <- function(simplex) paste(sort(simplex), collapse = "-")

#' @keywords internal
#' Example:
#'   .zz_simplex_keys(list(c(1, 2), c(3, 2)))
#'   #> [1] "1-2" "2-3"
.zz_simplex_keys <- function(simplex_list) vapply(simplex_list, .zz_simplex_key, character(1))

#' @keywords internal
#' Example:
#'   .zz_add_chains(c("1-2" = 1, "1-3" = -1), c("1-3" = 1, "2-3" = 1))
#'   #> 1-2 2-3
#'   #>   1   1
#'   # (the "1-3" entries cancel out and are dropped from the result)
.zz_add_chains <- function(a, b) {

  keys <- union(names(a), names(b))
  out <- setNames(numeric(length(keys)), keys)

  if (length(a) > 0) out[names(a)] <- out[names(a)] + a
  if (length(b) > 0) out[names(b)] <- out[names(b)] + b
  out[abs(out) > 1e-8]
}

#' @keywords internal
#' Example:
#'   .zz_chain_keys_union(list(c("1-2" = 1, "1-3" = -1), c("2-3" = 1)))
#'   #> [1] "1-2" "1-3" "2-3"
.zz_chain_keys_union <- function(vecs) unique(unlist(lapply(vecs, names)))

#' Full closure (every face of every dimension) of a maximal-simplex list
#' @keywords internal
#' Example:
#'   .zz_simplicial_closure(list(c(1, 2, 3)))
#'   #> list(1, 2, 3, c(1,2), c(1,3), c(2,3), c(1,2,3))
#'   # (every vertex, every edge, and the triangle itself)
.zz_simplicial_closure <- function(complex) {

  if (length(complex) == 0) return(list())
  max_dim <- max(sapply(complex, length)) - 1
  unlist(lapply(0:max_dim, function(p) faces(complex, p)), recursive = FALSE)
}

#' @keywords internal
#' Example:
#'   .zz_chains_to_matrix(list(c("1-2" = 1, "1-3" = -1), c("2-3" = 1)), keys = c("1-2", "1-3", "2-3"))
#'   #>     [,1] [,2]
#'   #> 1-2    1    0
#'   #> 1-3   -1    0
#'   #> 2-3    0    1
.zz_chains_to_matrix <- function(vecs, keys) {

  M <- matrix(0, nrow = length(keys), ncol = length(vecs))
  rownames(M) <- keys

  for (j in seq_along(vecs)) {
    v <- vecs[[j]]
    if (length(v) > 0) M[names(v), j] <- as.numeric(v)
  }
  M
}

#' Row-reduce (partial pivoting) a matrix, used internally for zigzag linear algebra
#' @keywords internal
#' Example (boundary of edges 1-2, 1-3, 2-3; rows = vertices 1, 2, 3):
#'   M <- matrix(c(1,1,0, -1,0,1, 0,-1,-1), nrow = 3, byrow = TRUE)
#'   .zz_gauss_jordan_eliminate(M)
#'   #> $R
#'   #>   [,1] [,2] [,3]
#'   #> 1    1    0   -1
#'   #> 2    0    1    1
#'   #> 3    0    0    0
#'   #> $pivots
#'   #> [1] 1 2
.zz_gauss_jordan_eliminate <- function(M, tol = 1e-8) {

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

#' Rank of a real matrix.
#' @keywords internal
#' Example:
#'   M <- matrix(c(1,1,0, -1,0,1, 0,-1,-1), nrow = 3, byrow = TRUE)
#'   .zz_matrix_rank(M)
#'   #> [1] 2
.zz_matrix_rank <- function(M) {

  if (is.null(M) || ncol(M) == 0 || nrow(M) == 0) return(0L)
  as.numeric(safe_rank(M))
}

#' Solve target = sum(c_i * gens[[i]]) over the reals; NULL if infeasible
#' @keywords internal
#' Example:
#'   target <- c("1-2" = 1, "2-3" = -1)
#'   gens <- list(c("1-2" = 1, "1-3" = -1), c("1-3" = 1, "2-3" = -1))
#'   .zz_solve_linear_combination(target, gens)
#'   #> [1] 1 1
#'   # (target = 1*gens[[1]] + 1*gens[[2]]; NULL would mean no such combination exists)
.zz_solve_linear_combination <- function(target, gens, tol = 1e-8) {

  keys <- union(names(target), .zz_chain_keys_union(gens))

  if (length(keys) == 0) return(if (length(gens) == 0) numeric(0) else rep(0, length(gens)))
  if (length(gens) == 0) {
    tvec <- .zz_chains_to_matrix(list(target), keys)[, 1]
    if (all(abs(tvec) < tol)) return(numeric(0)) else return(NULL)
  }
  Gmat <- .zz_chains_to_matrix(gens, keys)
  tvec <- .zz_chains_to_matrix(list(target), keys)[, 1]
  aug <- cbind(Gmat, tvec)
  rr <- .zz_gauss_jordan_eliminate(aug, tol)
  R <- rr$R; piv <- rr$pivots
  last_col <- ncol(aug)

  if (last_col %in% piv) return(NULL)
  coeffs <- numeric(ncol(Gmat))
  for (r in seq_along(piv)) {
    pc <- piv[r]
    if (pc <= ncol(Gmat)) coeffs[pc] <- R[r, last_col]
  }
  coeffs
}

#' Boundary of a single simplex as a named (by facet key) vector
#' @keywords internal
#' Example:
#'   .zz_simplex_boundary(c(1, 2, 3))
#'   #> 1-2 1-3 2-3
#'   #>   1  -1   1
.zz_simplex_boundary <- function(sx) {

  q <- length(sx) - 1
  if (q <= 0) return(setNames(numeric(0), character(0)))
  Mat <- suppressMessages(boundary(list(sx), q))
  row_faces <- faces(list(sx), q - 1)
  row_keys <- .zz_simplex_keys(row_faces)
  v <- as.numeric(as.matrix(Mat)[, 1])
  names(v) <- row_keys
  v[abs(v) > 1e-8]
}

#' Basis for the column space of a boundary matrix, as a list of named vectors
#' @keywords internal
#' Example (the single 2-simplex 1-2-3, boundary onto its 3 edges):
#'   Mat <- boundary(list(c(1, 2, 3)), 2)
#'   .zz_boundaries_basis(Mat, row_keys = c("1-2", "1-3", "2-3"))
#'   #> [[1]]
#'   #> 1-2 1-3 2-3
#'   #>   1  -1   1
.zz_boundaries_basis <- function(Mat, row_keys) {

  if (is.null(Mat) || ncol(Mat) == 0 || nrow(Mat) == 0) return(list())
  Md <- as.matrix(Mat)
  rr <- .zz_gauss_jordan_eliminate(Md)
  lapply(rr$pivots, function(pc) {
    v <- Md[, pc]; names(v) <- row_keys
    v[abs(v) > 1e-8]
  })
}

#' Basis for the null space of a boundary matrix (columns indexed by col_keys)
#' @keywords internal
#' Example (edges 1-2, 1-3, 2-3, boundary onto vertices 1, 2, 3):
#'   Mat <- boundary(list(c(1, 2), c(1, 3), c(2, 3)), 1)
#'   .zz_cycles_basis(Mat, col_keys = c("1-2", "1-3", "2-3"))
#'   #> [[1]]
#'   #> 1-2 1-3 2-3
#'   #>   1  -1   1
#'   # (the hollow triangle's single 1-cycle)
.zz_cycles_basis <- function(Mat, col_keys) {

  nc <- length(col_keys)
  if (nc == 0) return(list())
  if (is.null(Mat) || nrow(as.matrix(Mat)) == 0) {
    # zero map: every basis column is already a cycle
    return(lapply(seq_len(nc), function(j) setNames(1, col_keys[j])))
  }
  Md <- as.matrix(Mat)
  rr <- .zz_gauss_jordan_eliminate(Md)
  piv <- rr$pivots; R <- rr$R
  free <- setdiff(seq_len(nc), piv)
  lapply(free, function(fcol) {
    v <- numeric(nc)
    v[fcol] <- 1
    for (r in seq_along(piv)) v[piv[r]] <- -R[r, fcol]
    names(v) <- col_keys
    v[abs(v) > 1e-8]
  })
}

#' Extend B_basis to a basis of span(Z_basis); returns the added (complement) vectors
#' @keywords internal
#' Example (continuing the cycle from .zz_cycles_basis, no boundaries yet):
#'   Zp <- list(c("1-2" = 1, "1-3" = -1, "2-3" = 1))
#'   .zz_homology_basis(Z_basis = Zp, B_basis = list())
#'   #> [[1]]
#'   #> 1-2 1-3 2-3
#'   #>   1  -1   1
#'   # (nothing to cancel it against, so the cycle itself represents H1)
.zz_homology_basis <- function(Z_basis, B_basis) {

  if (length(Z_basis) == 0) return(list())
  keys <- unique(c(.zz_chain_keys_union(B_basis), .zz_chain_keys_union(Z_basis)))
  current <- if (length(B_basis) > 0) .zz_chains_to_matrix(B_basis, keys) else matrix(0, nrow = length(keys), ncol = 0)
  r_now <- .zz_matrix_rank(current)
  complement <- list()

  for (z in Z_basis) {
    zc <- .zz_chains_to_matrix(list(z), keys)
    test <- cbind(current, zc)
    r_new <- .zz_matrix_rank(test)
    if (r_new > r_now) {
      complement[[length(complement) + 1]] <- z
      current <- test
      r_now <- r_new
    }
  }
  complement
}

#' @keywords internal
#' Example:
#'   .zz_make_generator(c("1-2" = 1, "1-3" = -1, "2-3" = 1), birth = 0, id = 1)
#'   #> $vec
#'   #> 1-2 1-3 2-3
#'   #>   1  -1   1
#'   #> $birth
#'   #> [1] 0
#'   #> $id
#'   #> [1] 1
.zz_make_generator <- function(vec, birth, id) list(vec = vec, birth = birth, id = id)

#' Initialize the representative-cycle basis (actives) for every dimension of a complex
#' @keywords internal
#' Example (hollow triangle K0 = edges 1-2, 1-3, 2-3; one component, one loop):
#'   .zz_init_generators(list(c(1, 2), c(1, 3), c(2, 3)), birth_index = 0)
#'   #> $`0`: one generator, vec = c("1" = 1), birth = 0 (H0 = 1 component)
#'   #> $`1`: one generator, vec = c("1-2"=1,"1-3"=-1,"2-3"=1), birth = 0  (H1 = 1 loop)
#'   #> attr(, "next_id"): 3
.zz_init_generators <- function(complex, birth_index) {

  actives <- list()
  if (length(complex) == 0) {
    attr(actives, "next_id") <- 1L
    return(actives)
  }

  max_dim <- max(sapply(complex, length)) - 1
  next_id <- 1

  for (p in 0:max_dim) {
    p_faces <- faces(complex, p)
    if (length(p_faces) == 0) next
    p_keys <- .zz_simplex_keys(p_faces)
    Zp <- if (p == 0) {
      lapply(p_keys, function(k) setNames(1, k))
    } else {
      Mp <- suppressMessages(boundary(complex, p))
      .zz_cycles_basis(Mp, p_keys)
    }
    Mp1 <- suppressMessages(boundary(complex, p + 1))
    Bp <- .zz_boundaries_basis(Mp1, p_keys)
    Hp_basis <- .zz_homology_basis(Zp, Bp)
    gens <- lapply(Hp_basis, function(v) {
      g <- .zz_make_generator(v, birth_index, next_id)
      next_id <<- next_id + 1
      g
    })
    if (length(gens) > 0) actives[[as.character(p)]] <- gens
  }
  attr(actives, "next_id") <- next_id
  actives
}

#' Insert a single new simplex sx into the tracked zigzag state; may record a death (dim q-1)
#' or a birth (dim q). Mutates and returns the state list.
#' @keywords internal
#' Example (inserting the 2-simplex 1-2-3 into the hollow triangle above fills
#' the hole, so the H1 generator from .zz_init_generators dies immediately):
#'   state <- list(complex = list(c(1,2), c(1,3), c(2,3)), actives = init,
#'                 bars = list(), next_id = 3, from = 0L, to = 1L)
#'   state2 <- .zz_insert_simplex(state, c(1, 2, 3))
#'   state2$bars
#'   #> [[1]] list(dim = 1, birth = 0, death = 0)
#'   state2$actives[["1"]]
#'   #> list()   (no H1 generator survives - the loop was just filled in)
.zz_insert_simplex <- function(state, sx) {

  q <- length(sx) - 1
  key_sx <- .zz_simplex_key(sx)

  old_q <- faces(state$complex, q)
  old_q_keys <- .zz_simplex_keys(old_q)
  old_q_boundaries <- if (q > 0) lapply(old_q, .zz_simplex_boundary) else list()

  state$complex <- c(state$complex, list(sx))

  if (q == 0) {
    gen <- .zz_make_generator(setNames(1, key_sx), state$to, state$next_id)
    state$next_id <- state$next_id + 1
    state$actives[["0"]] <- c(state$actives[["0"]], list(gen))
    return(state)
  }

  dsx <- .zz_simplex_boundary(sx)
  sol <- .zz_solve_linear_combination(dsx, old_q_boundaries)

  if (!is.null(sol)) {
    # sx's boundary is already achievable from existing q-simplices -> a new q-cycle is born
    gen_vec <- setNames(1, key_sx)
    nz <- which(abs(sol) > 1e-8)
    if (length(nz) > 0) gen_vec <- .zz_add_chains(gen_vec, setNames(-sol[nz], old_q_keys[nz]))
    gen <- .zz_make_generator(gen_vec, state$to, state$next_id)
    state$next_id <- state$next_id + 1
    key_p <- as.character(q)
    state$actives[[key_p]] <- c(state$actives[[key_p]], list(gen))
  } else {
    # dsx is a genuinely new element of B_{q-1}: it must kill an existing (q-1) class
    key_pm1 <- as.character(q - 1)
    actives_qm1 <- state$actives[[key_pm1]]
    gens_for_solve <- c(old_q_boundaries, lapply(actives_qm1, `[[`, "vec"))
    sol2 <- .zz_solve_linear_combination(dsx, gens_for_solve)
    if (is.null(sol2)) stop("zigzag internal error: boundary of inserted simplex not in Z_{q-1}; invariant broken.")
    m <- length(old_q_boundaries)
    active_coeffs <- if (length(actives_qm1) > 0) sol2[(m + 1):length(sol2)] else numeric(0)
    involved <- which(abs(active_coeffs) > 1e-8)
    if (length(involved) == 0) stop("zigzag internal error: expected an active generator to die but found none.")
    ord <- involved[order(-sapply(actives_qm1[involved], `[[`, "birth"),
                           -sapply(actives_qm1[involved], `[[`, "id"))]
    lambda <- ord[1]
    died <- actives_qm1[[lambda]]
    state$bars[[length(state$bars) + 1]] <- list(dim = q - 1, birth = died$birth, death = state$from)
    state$actives[[key_pm1]] <- actives_qm1[-lambda]
  }
  state
}

#' Delete a single (currently maximal) simplex sx from the tracked zigzag state; may record
#' a death (dim q) or a birth (dim q-1). Mutates and returns the state list.
#' @keywords internal
#' Example (deleting edge 2-3 from the hollow triangle breaks the loop into a
#' path, so the H1 generator dies too - just via a different mechanism):
#'   state <- list(complex = list(c(1,2), c(1,3), c(2,3)), actives = init,
#'                 bars = list(), next_id = 3, from = 1L, to = 2L)
#'   state4 <- .zz_delete_simplex(state, c(2, 3))
#'   state4$bars
#'   #> [[1]] list(dim = 1, birth = 0, death = 1)
#'   state4$complex
#'   #> list(c(1,2), c(1,3))   (edge 2-3 removed)
.zz_delete_simplex <- function(state, sx) {

  q <- length(sx) - 1
  key_sx <- .zz_simplex_key(sx)

  # dim q: does removing sx force an active q-generator to disappear?
  key_p <- as.character(q)
  actives_q <- state$actives[[key_p]]
  if (length(actives_q) > 0) {
    coeffs <- sapply(actives_q, function(g) {
      v <- g$vec
      if (key_sx %in% names(v)) v[[key_sx]] else 0
    })
    S <- which(abs(coeffs) > 1e-8)
    if (length(S) > 0) {
      ord <- S[order(-sapply(actives_q[S], `[[`, "birth"), -sapply(actives_q[S], `[[`, "id"))]
      lambda <- ord[1]
      died <- actives_q[[lambda]]
      state$bars[[length(state$bars) + 1]] <- list(dim = q, birth = died$birth, death = state$from)
      for (k in S) {
        if (k == lambda) next
        ratio <- coeffs[k] / coeffs[lambda]
        actives_q[[k]]$vec <- .zz_add_chains(actives_q[[k]]$vec, -ratio * died$vec)
      }
      actives_q <- actives_q[-lambda]
    }
    state$actives[[key_p]] <- actives_q
  }

  # dim q-1: was sx essential to B_{q-1} (i.e. does removing it shrink the boundary space)?
  if (q >= 1) {
    other_q <- Filter(function(s) .zz_simplex_key(s) != key_sx, faces(state$complex, q))
    dsx <- .zz_simplex_boundary(sx)
    other_boundaries <- lapply(other_q, .zz_simplex_boundary)
    sol <- .zz_solve_linear_combination(dsx, other_boundaries)
    if (is.null(sol)) {
      # sx's boundary was independent: removing it frees a new (q-1) class
      key_pm1 <- as.character(q - 1)
      gen <- .zz_make_generator(dsx, state$to, state$next_id)
      state$next_id <- state$next_id + 1
      state$actives[[key_pm1]] <- c(state$actives[[key_pm1]], list(gen))
    }
  }

  state$complex <- Filter(function(s) .zz_simplex_key(s) != key_sx, state$complex)
  state
}

#' Compute the persistence barcode of a zigzag filtration of simplicial complexes
#'
#' Internally, each \eqn{K_i \to K_{i+1}} step is expanded into a sequence of single simplex insertions (faces before cofaces)
#' or single simplex deletions (cofaces before faces), exactly as \code{\link{boundary_info}}/\code{\link{persistence_pairs}}
#' assume for ordinary filtrations. At every elementary insertion or deletion, the representative-cycle basis of each affected
#' homology dimension is updated directly via linear algebra (reusing \code{\link{boundary}} and \code{\link{faces}} for every
#' boundary-matrix computation), an insertion either creates a new cycle (birth in dimension \eqn{q}) or turns an existing
#' cycle into a boundary (death in dimension \eqn{q-1}); a deletion either destroys an existing cycle (death in dimension \eqn{q})
#' or frees a previously-trivial cycle from being a boundary (birth in dimension \eqn{q-1}), where \eqn{q} is the
#' dimension of the simplex being inserted/deleted.
#'
#' @param complexes A list of length \eqn{n+1}: \code{complexes[[i+1]]} is \eqn{K_i}, a
#'   list of simplices exactly like the \code{simplices} argument of
#'   \code{\link{boundary}}/\code{\link{faces}} (each simplex a numeric vector; only the
#'   maximal simplices need to be listed, faces are inferred).
#' @param max_dimension Optional integer cap: dimensions above this are not tracked
#'   (saves work for large complexes where only e.g. \eqn{H_0}/\eqn{H_1} are of interest).
#'   \code{NULL} (default) tracks every dimension present.
#' @return A data frame with columns \code{dim}, \code{birth}, \code{death} (integer
#'   indices into \code{0, ..., n}, i.e. into \code{complexes}), one row per bar. Every
#'   bar uses a \strong{closed interval}: \code{death} is the last index at which
#'   the class is still present, and a class still alive at \eqn{K_n} is reported with
#'   \code{death = n}, since a zigzag filtration, unlike an ordinary one, has no
#'   canonical "infinity" to extend to. This differs by one from
#'   \code{\link{persistence_pairs}}'s convention, where \code{death} is the index of the
#'   killing simplex and the interval is half-open (\code{birth <= i < death});
#'   on a purely-growing filtration the two agree after \code{death_here = death_there - 1}.
#'
#' @export
#' @examples
#' # The K0..K4 triangle example: a hollow triangle whose 2-face is filled in and
#' # removed again, then broken by deleting an edge.
#' tri_edges <- list(c(1, 2), c(1, 3), c(2, 3))
#' K0 <- tri_edges
#' K1 <- c(tri_edges, list(c(1, 2, 3)))
#' K2 <- tri_edges
#' K3 <- c(tri_edges, list(4))
#' K4 <- list(c(1, 2), c(1, 3), 4)
#' bars <- zigzag_persistence(list(K0, K1, K2, K3, K4))
#' bars[bars$dim == 1, ]
zigzag_persistence <- function(complexes, max_dimension = NULL) {

  if (length(complexes) < 2) stop("complexes must have at least two complexes (K0 and K1).")
  n <- length(complexes) - 1

  init <- .zz_init_generators(complexes[[1]], birth_index = 0)
  state <- list(
    complex = complexes[[1]],
    actives = init,
    bars = list(),
    next_id = attr(init, "next_id"),
    from = 0L,
    to = 0L
  )
  # K0 -> K1 -> K2 -> ... -> Kn
  for (s in seq_len(n)) {
    closA <- .zz_simplicial_closure(complexes[[s]])
    closB <- .zz_simplicial_closure(complexes[[s + 1]])
    kA <- .zz_simplex_keys(closA)
    kB <- .zz_simplex_keys(closB)

    added <- closB[!(kB %in% kA)]
    removed <- closA[!(kA %in% kB)]

    if (length(added) > 0 && length(removed) > 0) {
      stop(sprintf("Step K%d -> K%d both adds and removes simplices; it is not a pure forward or backward inclusion. Split it into two steps.", s - 1, s))
    }

    state$from <- s - 1L # death
    state$to <- s # birth

    if (length(added) > 0) {
      ord <- order(sapply(added, length))
      for (sx in added[ord]) state <- .zz_insert_simplex(state, sx)
    } else if (length(removed) > 0) {
      ord <- order(-sapply(removed, length))
      for (sx in removed[ord]) state <- .zz_delete_simplex(state, sx)
    }
  }

  for (p_key in names(state$actives)) {
    for (g in state$actives[[p_key]]) {
      state$bars[[length(state$bars) + 1]] <- list(dim = as.integer(p_key), birth = g$birth, death = n)
    }
  }

  if (length(state$bars) == 0) {
    return(data.frame(dim = integer(0), birth = integer(0), death = integer(0)))
  }

  df <- data.frame(
    dim = as.integer(vapply(state$bars, `[[`, 0, "dim")),
    birth = as.integer(vapply(state$bars, `[[`, 0, "birth")),
    death = as.integer(vapply(state$bars, `[[`, 0, "death"))
  )
  if (!is.null(max_dimension)) df <- df[df$dim <= max_dimension, , drop = FALSE]
  df <- df[order(df$dim, df$birth, df$death), ]
  df
}
