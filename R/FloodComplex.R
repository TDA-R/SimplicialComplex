# ---------------------------------------------------------------------------
# Flood Complex for the SimplicialComplex package
#
# Pure-R re-implementation of the "flooder" Python package:
#   Graf, Pellizzoni, Uray, Huber, Kwitt (NeurIPS 2025),
#   "The Flood Complex: Large-Scale Persistent Homology on Millions of Points"
#   https://github.com/plus-rkwitt/flooder  (MIT)
#
# Pipeline:
#   1. Select landmarks via Farthest-Point Sampling (FPS).
#   2. Build the Delaunay triangulation on the landmarks (geometry::delaunayn).
#   3. For every top-dimensional Delaunay simplex, place a barycentric grid of
#      sample points on it. The filtration value ("flood time") of a simplex is
#      the covering radius: the maximum, over its grid points, of the distance
#      to the nearest witness point. Faces inherit their values from the grid
#      points that lie on them.
#   4. Make the filtration non-decreasing, then compute persistent homology.
#
# Dependencies:
#   geometry (Delaunay), RANN (kd-tree NN queries)
#   Optional: torch (GPU acceleration of the distance computation)
# ---------------------------------------------------------------------------

#' Farthest-Point Sampling of landmarks
#'
#' Selects \code{n_lms} landmarks from a point cloud via (exact) Farthest-Point
#' Sampling. Equivalent to \code{flooder::generate_landmarks} (which uses an
#' approximate bucket-FPS; the exact version below gives the same qualitative
#' coverage).
#'
#' @param points A numeric matrix (N x d) point cloud.
#' @param n_lms Number of landmarks to sample (<= N).
#' @param start_idx Index of the starting point. Defaults to 1
#'   (flooder defaults to index 0, i.e. the same first point).
#' @return A list with \code{landmarks} (n_lms x d matrix) and
#'   \code{indices} (row indices into \code{points}).
#'
#' @export
#' @examples
#' \dontrun{
#' pts <- matrix(rnorm(2000), ncol = 2)
#' lms <- generate_landmarks(pts, 50)
#' }
generate_landmarks <- function(points, n_lms, start_idx = 1) {
  points <- as.matrix(points)
  n <- nrow(points)
  if (n_lms <= 0) stop("Number of landmarks must be positive")
  n_lms <- min(n_lms, n)

  idx <- integer(n_lms)
  idx[1] <- start_idx
  # squared distance of every point to the current landmark set
  d2 <- rowSums(sweep(points, 2, points[start_idx, ])^2)
  if (n_lms > 1) {
    for (i in 2:n_lms) {
      idx[i] <- which.max(d2)
      d2_new <- rowSums(sweep(points, 2, points[idx[i], ])^2)
      d2 <- pmin(d2, d2_new)
    }
  }
  list(landmarks = points[idx, , drop = FALSE], indices = idx)
}

#' @keywords internal
#' Barycentric grid on the unit simplex (replicates flooder's generate_grid).
#' Returns the weight matrix (C x (dim+1)) plus, for every face of the
#' standard simplex, which local vertices span it and which grid rows lie on it.
.flood_grid <- function(n, dim) {
  # all non-negative integer vectors of length dim+1 summing to n-1
  if (dim == 0) {
    grid <- matrix(n - 1, 1, 1)
  } else {
    combs <- t(utils::combn(0:(n + dim - 2), dim)) # C x dim
    padded <- cbind(-1, combs, n + dim - 1) # C x (dim+2)
    grid <- matrix(t(apply(padded, 1, diff)) - 1, ncol = dim + 1) # C x (dim+1)
  }
  weights <- grid / (n - 1)

  face_rows <- list() # rows of the grid lying on each face
  face_vertices <- list() # local vertex indices (1-based) spanning each face
  m <- 0
  for (k in 0:dim) { # k = number of barycentric coordinates forced to zero
    zero_sets <- utils::combn(seq_len(dim + 1), k, simplify = FALSE)
    if (k == 0) zero_sets <- list(integer(0))
    for (zs in zero_sets) {
      m <- m + 1
      if (length(zs) == 0) {
        rows <- seq_len(nrow(grid))
      } else {
        rows <- which(rowSums(grid[, zs, drop = FALSE] != 0) == 0)
      }
      face_rows[[m]] <- rows
      face_vertices[[m]] <- setdiff(seq_len(dim + 1), zs)
    }
  }
  list(weights = weights, face_rows = face_rows, face_vertices = face_vertices)
}

#' @keywords internal
#' Nearest-neighbour distances, CPU backend (kd-tree via RANN, with a
#' base-R brute-force fallback when RANN is not installed).
.flood_nn_cpu <- function(points, queries) {
  if (requireNamespace("RANN", quietly = TRUE)) {
    RANN::nn2(data = points, query = queries, k = 1)$nn.dists[, 1]
  } else {
    .flood_nn_base(points, queries)
  }
}

#' @keywords internal
#' Chunked brute-force nearest-neighbour distances in base R.
.flood_nn_base <- function(points, queries, chunk = 2048L) {
  pt2 <- rowSums(points^2)
  nq <- nrow(queries)
  out <- numeric(nq)
  for (s in seq(1L, nq, by = chunk)) {
    e <- min(nq, s + chunk - 1L)
    Q <- queries[s:e, , drop = FALSE]
    d2 <- outer(rowSums(Q^2), pt2, "+") - 2 * (Q %*% t(points))
    j <- max.col(-d2, ties.method = "first") # fast row-wise argmin
    out[s:e] <- sqrt(pmax(d2[cbind(seq_len(nrow(d2)), j)], 0))
  }
  out
}

#' @keywords internal
#' Nearest-neighbour distances, torch backend (GPU if CUDA is available).
#' Brute-force cdist in chunks; mirrors flooder's non-Triton CUDA path.
.flood_nn_torch <- function(points_t, queries, device, chunk = 4096L) {
  q <- torch::torch_tensor(queries, dtype = torch::torch_float32(), device = device)
  nq <- nrow(queries)
  out <- numeric(nq)
  for (s in seq(1L, nq, by = chunk)) {
    e <- min(nq, s + chunk - 1L)
    idx <- torch::torch_tensor(as.integer(s:e), dtype = torch::torch_long(),
                               device = device)
    d <- torch::torch_cdist(q$index_select(1L, idx), points_t)
    out[s:e] <- as.numeric(torch::torch_min(d, dim = 2L)[[1]]$cpu())
  }
  out
}

#' Construct a Flood complex
#'
#' Builds the Flood complex of a point cloud: a Delaunay complex on a set of
#' landmarks whose simplices are filtered by their covering radius with
#' respect to the full point cloud ("flood time").
#'
#' @param points A numeric matrix (N x d) of witness points.
#' @param landmarks Either an integer (number of FPS landmarks) or a numeric matrix (N_l x d) of explicit landmark coordinates.
#' @param max_dimension Top dimension of the simplices. Defaults to d.
#' @param points_per_edge Grid resolution per simplex edge (accuracy vs. speed trade-off). Defaults to 30, as in flooder.
#' @param backend One of \code{"auto"}, \code{"cpu"}, \code{"torch"}.
#'   \code{"cpu"} uses a kd-tree (RANN). \code{"torch"} uses the torch R package and runs on CUDA when available. \code{"auto"} picks torch only if a CUDA device is present, else the kd-tree.
#' @param batch_points Maximum number of grid points processed per batch (bounds memory). Defaults to 2^22.
#' @param delaunay Optional precomputed Delaunay triangulation of the landmarks: an m x (d+1) integer matrix of 1-based landmark indices. If NULL (default), computed via \code{geometry::delaunayn}.
#' @return A list of class \code{"flood_complex"} with elements
#'   \code{simplices} (list of integer vectors, landmark indices),
#'   \code{filtration} (numeric vector of flood times),
#'   \code{landmarks} (matrix), \code{landmark_indices} (or NULL).
#'
#' @export
#' @examples
#' \dontrun{
#' pts <- matrix(rnorm(3000), ncol = 2)
#' fc <- flood_complex(pts, landmarks = 40)
#' }
flood_complex <- function(points,
                          landmarks,
                          max_dimension = NULL,
                          points_per_edge = 30,
                          backend = c("auto", "cpu", "torch"),
                          batch_points = 2^22,
                          delaunay = NULL) {
  backend <- match.arg(backend)
  points <- as.matrix(points)
  d_amb <- ncol(points)
  if (is.null(max_dimension)) max_dimension <- d_amb

  lm_idx <- NULL
  if (length(landmarks) == 1 && is.numeric(landmarks)) {
    lm <- generate_landmarks(points, as.integer(landmarks))
    landmarks <- lm$landmarks
    lm_idx <- lm$indices
  } else {
    landmarks <- as.matrix(landmarks)
  }
  if (ncol(landmarks) != d_amb) stop("landmarks and points must have the same dimension")

  # backend selection
  use_torch <- FALSE
  device <- NULL
  if (backend == "torch" ||
      (backend == "auto" &&
       requireNamespace("torch", quietly = TRUE) &&
       torch::cuda_is_available())) {
    if (!requireNamespace("torch", quietly = TRUE)) {
      stop("backend = 'torch' requested but the torch package is not installed")
    }
    use_torch <- TRUE
    device <- if (torch::cuda_is_available()) "cuda" else "cpu"
    points_t <- torch::torch_tensor(points, dtype = torch::torch_float32(),
                                    device = device)
  }

  # Delaunay triangulation on the landmarks
  if (is.null(delaunay)) {
    if (!requireNamespace("geometry", quietly = TRUE)) {
      stop("Package 'geometry' is required (or pass a precomputed 'delaunay')")
    }
    delaunay <- geometry::delaunayn(landmarks, options = "Qt Qbb Qc Qz")
  }
  top <- matrix(t(apply(delaunay, 1, sort)), ncol = ncol(delaunay)) # m x (d_amb+1)

  # If max_dimension < ambient dimension, use the unique max_dimension-faces
  # of the Delaunay top simplices as top simplices.
  if (max_dimension < d_amb) {
    sub <- utils::combn(seq_len(d_amb + 1), max_dimension + 1, simplify = FALSE)
    top <- unique(do.call(rbind, lapply(sub, function(s) top[, s, drop = FALSE])))
  }
  d_top <- ncol(top) - 1
  m <- nrow(top)

  # grid on the standard simplex
  g <- .flood_grid(points_per_edge, d_top)
  W <- g$weights # C x (d_top+1)
  C <- nrow(W)

  # covering radii
  filt <- new.env(hash = TRUE, parent = emptyenv())
  batch_size <- max(1L, as.integer(batch_points / C))

  for (s0 in seq(1L, m, by = batch_size)) {
    s1 <- min(m, s0 + batch_size - 1L)
    nb <- s1 - s0 + 1L
    Sb <- top[s0:s1, , drop = FALSE]

    # grid points on each simplex of the batch: (nb*C) x d_amb
    P <- matrix(0, nb * C, d_amb)
    for (i in seq_len(nb)) {
      P[((i - 1L) * C + 1L):(i * C), ] <- W %*% landmarks[Sb[i, ], , drop = FALSE]
    }

    dists <- if (use_torch) .flood_nn_torch(points_t, P, device)
             else .flood_nn_cpu(points, P)
    D <- matrix(dists, nrow = nb, ncol = C, byrow = TRUE)

    # covering radius of every face of every simplex in the batch
    for (f in seq_along(g$face_rows)) {
      rows <- g$face_rows[[f]]
      verts <- g$face_vertices[[f]]
      vals <- if (length(rows) == 1L) D[, rows] else {
        M <- D[, rows, drop = FALSE]
        M[cbind(seq_len(nrow(M)), max.col(M, ties.method = "first"))]
      }
      faces_b <- Sb[, verts, drop = FALSE]                  # nb x |face|
      keys <- apply(faces_b, 1, paste, collapse = " ")
      for (i in seq_len(nb)) assign(keys[i], vals[i], envir = filt)
    }
  }

  # collect + make filtration non-decreasing
  keys <- ls(filt)
  simplices <- lapply(strsplit(keys, " ", fixed = TRUE), as.integer)
  values <- unlist(mget(keys, envir = filt), use.names = FALSE)
  dims <- lengths(simplices) - 1L

  ord <- order(dims)
  simplices <- simplices[ord]; values <- values[ord]; dims <- dims[ord]
  key_of <- vapply(simplices, paste, "", collapse = " ")
  lookup <- new.env(hash = TRUE, parent = emptyenv())
  for (i in seq_along(key_of)) assign(key_of[i], i, envir = lookup)

  for (i in seq_along(simplices)) {
    if (dims[i] == 0L) next
    s <- simplices[[i]]
    facets <- utils::combn(s, length(s) - 1L)
    for (j in seq_len(ncol(facets))) {
      fi <- get(paste(facets[, j], collapse = " "), envir = lookup)
      if (values[fi] > values[i]) values[i] <- values[fi]
    }
  }

  structure(
    list(simplices = simplices, filtration = values, landmarks = landmarks, landmark_indices = lm_idx),
    class = "flood_complex"
  )
}

#' Flood filtration in SimplicialComplex format
#'
#' Wraps \code{\link{flood_complex}} and returns the filtration in the same
#' format as \code{\link{build_filtration}}: a list of
#' \code{list(simplex = <integer vector>, t = <numeric>)}, sorted by
#' (time, dimension, lexicographic order). The result plugs directly into
#' \code{boundary_info()} / \code{extract_persistence_pairs()} /
#' \code{plot_persistence()}, as well as into the faster
#' \code{\link{flood_persistence}}.
#'
#' @param points A numeric matrix (N x d) point cloud.
#' @param landmarks Number of FPS landmarks, or an explicit landmark matrix.
#' @param ... Passed on to \code{\link{flood_complex}}.
#' @return A filtration list compatible with \code{boundary_info}.
#'
#' @export
#' @examples
#' \dontrun{
#' pts <- matrix(rnorm(2000), ncol = 2)
#' filtration <- build_flood_filtration(pts, landmarks = 30)
#' pairs <- flood_persistence(filtration)
#' }
build_flood_filtration <- function(points, landmarks, ...) {
  as_filtration(flood_complex(points, landmarks, ...))
}

#' Convert a flood_complex object to a filtration list
#'
#' @param fc A \code{"flood_complex"} object.
#' @return A filtration list compatible with \code{boundary_info}.
#'
#' @export
as_filtration <- function(fc) {
  filist <- Map(function(s, t) list(simplex = s, t = t),
                fc$simplices, fc$filtration)
  ord <- order(vapply(filist, `[[`, 0, "t"),
               lengths(lapply(filist, `[[`, "simplex")),
               vapply(filist, function(x) paste(x$simplex, collapse = "-"), ""))
  filist[ord]
}

#' Persistence pairs via sparse boundary reduction
#'
#' Computes persistence pairs from a filtration list with the standard
#' column-reduction algorithm, but on sparse columns (integer index vectors
#' over GF(2)) instead of a dense matrix. Produces the same output format as
#' \code{extract_persistence_pairs(filist, res$last_1, res$pivot_owner)} while
#' scaling to the much larger complexes produced by
#' \code{\link{build_flood_filtration}}.
#'
#' @param filist A filtration list (from \code{build_flood_filtration} or
#'   \code{build_filtration}).
#' @param max_dimension Optional maximum homology dimension to report. When
#'   set, one extra dimension is kept internally so dimension-\code{max_dimension}
#'   classes still get correct death times from their true killers, and only
#'   that extra dimension is dropped from the output - see
#'   \code{\link{restrict_filtration}}'s Details, and
#'   \code{\link{persistence_pairs}} which applies the same correction.
#'   Leave \code{NULL} (default) to report every dimension present.
#' @return A data frame with columns \code{dim}, \code{birth}, \code{death}.
#'
#' @export
flood_persistence <- function(filist, max_dimension = NULL) {
  if (!is.null(max_dimension)) {
    filist <- restrict_filtration(filist, max_dimension + 1)
  }
  n <- length(filist)
  keys <- vapply(filist, function(x) paste(x$simplex, collapse = " "), "")
  index <- new.env(hash = TRUE, parent = emptyenv())
  for (i in seq_len(n)) assign(keys[i], i, envir = index)

  cols <- vector("list", n)
  for (i in seq_len(n)) {
    s <- filist[[i]]$simplex
    if (length(s) == 1L) { cols[[i]] <- integer(0); next }
    fmat <- utils::combn(s, length(s) - 1L)
    cols[[i]] <- sort(vapply(seq_len(ncol(fmat)), function(j)
      get(paste(fmat[, j], collapse = " "), envir = index), 0L))
  }

  symdiff <- function(a, b) sort.int(c(a[!(a %in% b)], b[!(b %in% a)]))

  pivot_owner <- rep(NA_integer_, n)
  for (j in seq_len(n)) {
    col <- cols[[j]]
    repeat {
      if (length(col) == 0L) break
      piv <- col[length(col)]
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
  # zero column (positive simplex) that is never paired
  zero_cols <- which(vapply(cols, length, 0L) == 0L)
  essential <- zero_cols[is.na(pivot_owner[zero_cols])]
  if (length(essential)) {
    pairs <- rbind(pairs, data.frame(
      dim = dims[essential], birth = ts[essential], death = Inf))
  }
  rownames(pairs) <- NULL
  pairs <- pairs[order(pairs$dim, pairs$birth), ]

  if (!is.null(max_dimension)) {
    pairs <- pairs[pairs$dim <= max_dimension, , drop = FALSE]
    rownames(pairs) <- NULL
  }
  pairs
}
