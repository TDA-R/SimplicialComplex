# Discrete Morse theory: persistence-guided discrete Morse vector fields (DMVF) and graph reconstruction from a density field.
#
# This follows Chapter 10 "Discrete Morse Theory and Applications" from "Computational Topology for Data Analysis":
#
#   - lower_star_filtration():  generalizes build_cubical_filtration() to any
#     triangulation, builds the simplex-wise lower-star filtration F_f of a
#     PL vertex function f.
#
#   - partial_pers_dmvf(): Algorithm 19 (SimplePersDMVF) combined with
#     the threshold-delta simplification: cancels
#     vertex-edge persistence pairs with persistence <= delta, producing a
#     discrete Morse vector field (DMVF) on a 1-complex (graph).
#
#   - morse_recon(): Algorithm 20 (MorseRecon): given a
#     triangulated 2-complex and a density function rho concentrated around
#     a hidden geometric graph, reconstructs the graph as the union of the
#     1-unstable manifolds ("mountain ridges") of f = -rho, after cancelling
#     low-persistence (noise) vertex-edge pairs. Internally refines the
#     persistence of edges that do NOT merge two components ("creator"
#     edges) using the edge-triangle (H1) pairing of the full 2-complex.
#
#   - collect_g(): Algorithm 21 (CollectG): collects the
#     1-unstable manifolds of the surviving high-persistence critical edges
#     into the reconstructed graph. Exposed separately so it can also be
#     used directly on a partial_pers_dmvf() result for a plain graph.

.uf_find <- function(uf, x) {
  root <- x
  while (uf$parent[root] != root) root <- uf$parent[root]
  while (uf$parent[x] != root) {
    nxt <- uf$parent[x]
    uf$parent[x] <- root
    x <- nxt
  }
  list(uf = uf, root = root)
}

.uf_union <- function(uf, rx, ry) {
  if (rx == ry) return(list(uf = uf, root = rx))
  if (uf$rank[rx] < uf$rank[ry]) { tmp <- rx; rx <- ry; ry <- tmp }
  uf$parent[ry] <- rx
  if (uf$rank[rx] == uf$rank[ry]) uf$rank[rx] <- uf$rank[rx] + 1L
  list(uf = uf, root = rx)
}

# lower-star filtration
#' Generalizes \code{\link{build_cubical_filtration}} to any triangulation:
#' given the maximal simplices of a simplicial complex and a function defined at its vertices, builds the simplex-wise lower-star filtration \eqn{\mathcal{F}_f}.
#'
#' @param top_simplices A list of maximal simplices (each an integer vector of 1-based vertex ids).
#' @param f A plain numeric vector giving the function value at each vertex; \code{f[v]} is the value at vertex \code{v},
#'   so vertex ids must be integers in \code{1:length(f)}.
#' @return A filtration list: one \code{list(simplex =integer vector of vertex ids, t = numeric)} per simplex (every face of every maximal simplex),
#'   sorted by \code{(t, dimension, lexicographic order)}.
#'
#' @importFrom utils combn
#' @export
#' @examples
#' # two triangles sharing an edge, function increasing away from vertex 1
#' triangles <- list(c(1, 2, 3), c(2, 3, 4))
#' f <- c(0, 1, 1, 2)
#' filtration <- lower_star_filtration(triangles, f)
#' length(filtration) # 4 vertices + 5 edges + 2 triangles
lower_star_filtration <- function(top_simplices, f) {

  stopifnot(is.numeric(f), length(top_simplices) > 0)

  seen <- new.env(hash = TRUE, parent = emptyenv())
  simplex_info <- list()

  add_simplex <- function(verts) {
    verts <- sort(unique(as.integer(verts)))
    key <- paste(verts, collapse = "-")
    if (exists(key, envir = seen, inherits = FALSE)) return(invisible(NULL))
    assign(key, TRUE, envir = seen)
    simplex_info[[length(simplex_info) + 1L]] <<- list(simplex = verts, t = max(f[verts]))
  }

  for (s in top_simplices) {
    d <- length(s)
    for (k in seq_len(d)) { # k is how many rows
      combs <- combn(s, k)
      for (j in seq_len(ncol(combs))) add_simplex(combs[, j])
    }
  }

  ord <- order(
    vapply(simplex_info, `[[`, 0, "t"),
    lengths(lapply(simplex_info, `[[`, "simplex")),
    vapply(simplex_info, function(x) paste(x$simplex, collapse = "-"), "")
  )
  simplex_info[ord]
}

# structural (H0) vertex-edge persistence pass
# Pure union-find pass over a 1-skeleton filtration: classifies every edge as a "destroyer" (merges two components; gets a finite persistence
# |f(e) - f(r)| where r is the younger of the two roots, the "elder rule") or a "creator" (both endpoints already in the same component; persistence
# Inf until possibly refined by a later H1 computation).
.vertex_edge_pers <- function(filist1) {

  is_v <- vapply(filist1, function(x) length(x$simplex) == 1L, logical(1))
  is_e <- vapply(filist1, function(x) length(x$simplex) == 2L, logical(1))

  if (!all(is_v | is_e)) {stop("filist1 must contain only vertices and edges, use restrict_filtration(filist, 1) first.")}

  vfil <- filist1[is_v]
  vids <- vapply(vfil, function(x) x$simplex[1], integer(1))
  n <- length(vids)

  idx <- stats::setNames(seq_len(n), as.character(vids)) # vertex id -> uf index
  arrival <- stats::setNames(seq_len(n), as.character(vids)) # filtration arrival order
  vt <- stats::setNames(vapply(vfil, `[[`, 0, "t"), as.character(vids))

  uf <- list(parent = seq_len(n), rank = integer(n))
  logical_root <- vids # logical_root[uf-root-index] = id of the earliest vertex in that component

  n_edges <- sum(is_e)
  out <- data.frame(
    u = integer(n_edges), v = integer(n_edges), t = numeric(n_edges),
    persistence = numeric(n_edges), type = character(n_edges),
    stringsAsFactors = FALSE
  )
  e_i <- 0L

  for (item in filist1) {
    s <- item$simplex
    if (length(s) != 2L) next
    u <- s[1]; v <- s[2]
    fu <- .uf_find(uf, idx[[as.character(u)]]); uf <- fu$uf; ru <- fu$root
    fv <- .uf_find(uf, idx[[as.character(v)]]); uf <- fv$uf; rv <- fv$root

    e_i <- e_i + 1L
    out$u[e_i] <- u; out$v[e_i] <- v; out$t[e_i] <- item$t

    if (ru == rv) {
      out$type[e_i] <- "creator"
      out$persistence[e_i] <- Inf
    } else {
      r1 <- logical_root[ru]; r2 <- logical_root[rv]
      if (arrival[[as.character(r1)]] < arrival[[as.character(r2)]]) {
        tmp <- r1; r1 <- r2; r2 <- tmp
      }
      # r1 = younger (later-arriving) root, r2 = elder root -> survives (elder rule)
      out$type[e_i] <- "destroyer"
      out$persistence[e_i] <- abs(item$t - vt[[as.character(r1)]])
      un <- .uf_union(uf, ru, rv); uf <- un$uf
      logical_root[un$root] <- r2
    }
  }
  out[seq_len(e_i), , drop = FALSE]
}

# edge-triangle (H1) persistence, for refining H0-creator edges
# Column-reduction restricted to the edge/triangle slice of a full 2-complex filtration.
.dim1_triangle_pairing <- function(filist2) {

  is_e <- vapply(filist2, function(x) length(x$simplex) == 2L, logical(1))
  is_t <- vapply(filist2, function(x) length(x$simplex) == 3L, logical(1))

  empty <- list(pairs = data.frame(u = integer(0), v = integer(0), persistence = numeric(0)))
  if (!any(is_t) || !any(is_e)) return(empty)

  red <- .reduce_gf2_boundary(filist2)
  pivot_owner <- red$pivot_owner

  killed <- which(is_e & !is.na(pivot_owner))
  paired <- killed[is_t[pivot_owner[killed]]]
  if (length(paired) == 0) return(empty)

  s <- do.call(rbind, lapply(filist2[paired], `[[`, "simplex"))
  list(pairs = data.frame(
    u = s[, 1], v = s[, 2],
    persistence = vapply(paired, function(row) {
      j <- pivot_owner[row]
      abs(filist2[[j]]$t - filist2[[row]]$t)
    }, numeric(1)
    )
  ))
}

# build the matching forest for a given threshold delta
#
# A tree/matching edge is exactly a "destroyer" edge with persistence <= delta: all "creator" edges,
# and all "destroyer" edges with persistence > delta, stay critical, this is exactly K^1 \ T.
# Roots are recomputed per surviving connected component as its earliest-arriving vertex;
# parent(v) is then read off a single BFS/traversal per component, matching v with the edge (v, parent(v))
.build_forest <- function(filist1, pers_df, delta) {

  is_v <- vapply(filist1, function(x) length(x$simplex) == 1L, logical(1))
  vfil <- filist1[is_v]
  vids <- vapply(vfil, function(x) x$simplex[1], integer(1))
  arrival <- stats::setNames(seq_along(vids), as.character(vids))

  is_tree <- pers_df$type == "destroyer" & pers_df$persistence <= delta
  tree <- pers_df[is_tree, , drop = FALSE]
  crit <- pers_df[!is_tree, , drop = FALSE]

  adj <- stats::setNames(vector("list", length(vids)), as.character(vids))

  for (i in seq_len(nrow(tree))) {
    uk <- as.character(tree$u[i]); vk <- as.character(tree$v[i])
    adj[[uk]] <- c(adj[[uk]], tree$v[i])
    adj[[vk]] <- c(adj[[vk]], tree$u[i])
  }

  parent <- stats::setNames(rep(NA_integer_, length(vids)), as.character(vids))
  visited <- stats::setNames(rep(FALSE, length(vids)), as.character(vids))
  roots <- integer(0)

  order_verts <- vids[order(arrival[as.character(vids)])]

  for (r in order_verts) {

    rk <- as.character(r)
    if (visited[[rk]]) next
    roots <- c(roots, r)
    queue <- r
    visited[[rk]] <- TRUE

    while (length(queue) > 0) {
      cur <- queue[1]; queue <- queue[-1]
      for (nb in adj[[as.character(cur)]]) {
        nk <- as.character(nb)
        if (!visited[[nk]]) {
          visited[[nk]] <- TRUE
          parent[[nk]] <- cur
          queue <- c(queue, nb)
        }
      }
    }
  }

  list(tree_edges = tree, critical_edges = crit, parent = parent, roots = roots, delta = delta)
}

# SimplePersDMVF

#' Persistence-guided discrete Morse vector field on a graph (1-complex)
#'
#' Implements Algorithm 19 (\code{SimplePersDMVF}) together
#' with the threshold-\eqn{\delta} simplification: vertex-edge persistence pairs with persistence
#' at most \code{delta} are cancelled, leaving a discrete Morse vector field
#' where every non-critical vertex \eqn{v} is matched with the tree edge
#' connecting it to its component's root, and every root is a critical vertex.
#'
#' @param filist1 A filtration list restricted to vertices and edges only (e.g. via \code{restrict_filtration(filist, 1)}),
#'    as produced by \code{\link{lower_star_filtration}} or \code{\link{build_filtration}}.
#' @param delta Persistence threshold; vertex-edge pairs with persistence \code{<= delta} are cancelled (become matched V-field arrows).
#'   Default \code{0} cancels only exactly-zero-persistence pairs; use a larger value to denoise more aggressively, or \code{Inf} to cancel
#'   every finite pair (maximal simplification, Theorem 10.5).
#' @return A list with components:
#'   \item{pers}{data.frame(u, v, t, persistence, type) for every edge, \code{type} is \code{"destroyer"} (merges two components) or \code{"creator"} (closes a cycle; persistence \code{Inf})}
#'   \item{tree_edges}{the subset of \code{pers} used as matching/DMVF edges}
#'   \item{critical_edges}{the rest of \code{pers}, i.e. \eqn{K^1 \setminus T} -- exactly the input \code{\link{collect_g}} needs}
#'   \item{parent}{named integer vector, \code{parent[["v"]]} is the vertex on \eqn{v}'s tree edge towards its component root (\code{NA} for roots)}
#'   \item{roots}{integer vector of critical (root) vertices}
#'   \item{delta}{the threshold used}
#' @export
partial_pers_dmvf <- function(filist1, delta = 0) {
  pers_df <- .vertex_edge_pers(filist1)
  forest <- .build_forest(filist1, pers_df, delta)
  c(list(pers = pers_df), forest)
}

# CollectG

#' Collect 1-unstable manifolds into a reconstructed graph
#'
#' Implements Algorithm 21 (\code{CollectG}): for every
#' critical edge \eqn{e = (u, v)} (i.e. every edge not used as a DMVF
#' matching/tree edge - by construction these all have persistence greater
#' than \code{dmvf$delta}, so unlike the book's pseudocode no extra filter
#' is needed here), unions \eqn{e} with the unique tree paths from \eqn{u}
#' and from \eqn{v} up to their respective roots.
#'
#' @param dmvf A result of \code{\link{partial_pers_dmvf}} (or the equivalent list built inside \code{\link{morse_recon}}).
#' @return A \code{data.frame(u, v)} of the (deduplicated, undirected) edges of the reconstructed graph.
#' @export
collect_g <- function(dmvf) {

  crit <- dmvf$critical_edges
  if (nrow(crit) == 0) return(data.frame(u = integer(0), v = integer(0)))

  path_to_root <- function(v) {
    edges <- list()
    cur <- v
    repeat {
      p <- dmvf$parent[[as.character(cur)]]
      if (is.na(p)) break
      edges[[length(edges) + 1L]] <- c(cur, p)
      cur <- p
    }
    edges
  }

  out <- vector("list", nrow(crit) * 3L) # rough upper bound per edge, grown as needed
  k <- 0L
  add <- function(pair) {k <<- k + 1L; out[[k]] <<- pair}

  for (i in seq_len(nrow(crit))) {
    add(c(crit$u[i], crit$v[i]))
    for (e in path_to_root(crit$u[i])) add(e)
    for (e in path_to_root(crit$v[i])) add(e)
  }
  out <- out[seq_len(k)]

  m <- do.call(rbind, out)
  m <- unique(t(apply(m, 1, sort)))
  data.frame(u = m[, 1], v = m[, 2])
}

# MorseRecon

#' Reconstruct a hidden graph from a density field
#'
#' Implements Algorithm 20 (\code{MorseRecon}):
#' given a triangulated domain and a density function \code{rho} that concentrates around a hidden geometric graph \eqn{G},
#' computes the "mountain ridges" of \eqn{f = -\rho} - the 1-unstable manifolds of the discrete gradient field after cancelling
#' vertex-edge persistence pairs with persistence at most \code{delta} - as an approximation \eqn{\hat G} of \eqn{G}.
#'
#' @param top_simplices List of maximal simplices (e.g. triangles) of the ambient 2-complex, as accepted by \code{\link{lower_star_filtration}}.
#' @param rho A plain numeric vector, the density value at each vertex.
#' @param delta Persistence threshold used to cancel low-persistence vertex-edge pairs (noise); larger values denoise more aggressively.
#' @param vertex_coords Optional n x 2 matrix of vertex coordinates, used only by \code{\link{plot_morse_recon}} for drawing.
#'
#' @return An object of class \code{"morse_recon"}: a list with
#'   \item{filtration}{the full lower-star filtration of the 2-complex}
#'   \item{dmvf}{the vertex-edge DMVF: \code{pers}, \code{tree_edges},
#'     \code{critical_edges}, \code{parent}, \code{roots}, \code{delta}
#'     (same shape as a \code{\link{partial_pers_dmvf}} result, but with
#'     H0-creator edges' persistence refined against triangles, see Details)}
#'   \item{graph_edges}{data.frame(u, v) - the edges of \eqn{\hat G}}
#'   \item{rho, delta, coords}{the inputs, kept for plotting/inspection}
#'
#' @importFrom stats setNames
#' @export
morse_recon <- function(top_simplices, rho, delta = 0, vertex_coords = NULL) {

  f <- -as.numeric(rho)

  filist2 <- lower_star_filtration(top_simplices, f)
  filist1 <- restrict_filtration(filist2, max_dimension = 1)

  pers_df <- .vertex_edge_pers(filist1)
  tri_pairing <- .dim1_triangle_pairing(filist2)

  if (nrow(tri_pairing$pairs) > 0) {

    creator_idx <- which(pers_df$type == "creator")

    if (length(creator_idx) > 0) {

      keys <- paste(pers_df$u[creator_idx], pers_df$v[creator_idx], sep = "-")
      pair_keys <- paste(tri_pairing$pairs$u, tri_pairing$pairs$v, sep = "-")

      m <- match(keys, pair_keys)
      hit <- !is.na(m)
      pers_df$persistence[creator_idx[hit]] <- tri_pairing$pairs$persistence[m[hit]]
      drop_idx <- creator_idx[hit][tri_pairing$pairs$persistence[m[hit]] <= delta]

      if (length(drop_idx) > 0) pers_df <- pers_df[-drop_idx, , drop = FALSE]
    }
  }

  forest <- .build_forest(filist1, pers_df, delta)
  dmvf <- c(list(pers = pers_df), forest)
  graph_edges <- collect_g(dmvf)

  structure(
    list(filtration = filist2, dmvf = dmvf, graph_edges = graph_edges, rho = rho, delta = delta, coords = vertex_coords),
    class = "morse_recon"
  )
}
