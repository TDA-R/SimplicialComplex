#' @keywords internal
diagram_points <- function(df, dimension) {
  sub <- df[df$dim == dimension, , drop = FALSE]
  cbind(birth = sub$birth, death = sub$death)
}

#' @keywords internal
cap_essential <- function(X, Y, dimension) {
  finite_vals <- c(X[is.finite(X)], Y[is.finite(Y)])
  cap <- if (length(finite_vals) > 0) 2 * max(finite_vals) else 1
  if (!is.finite(cap) || cap <= 0) cap <- 1

  n_essential <- sum(!is.finite(X)) + sum(!is.finite(Y))
  if (n_essential > 0) {
    message(sprintf(
      "Capping %d essential (death = Inf) pair(s) in dimension %d at death = %.4g so they are compared by birth time rather than discarded.",
      n_essential, dimension, cap
    ))
  }

  X[!is.finite(X)] <- cap
  Y[!is.finite(Y)] <- cap
  list(X = X, Y = Y)
}

#' @keywords internal
ground_dist <- function(X, Y, ground = c("L2", "Linf")) {
  ground <- match.arg(ground)
  if (ground == "L2") return(pairwise_dist(X, Y))
  # Linf (Chebyshev): max coordinate-wise difference
  n <- nrow(X); m <- nrow(Y)
  out <- matrix(0, n, m)
  for (k in seq_len(ncol(X))) out <- pmax(out, abs(outer(X[, k], Y[, k], "-")))
  out
}

#' @keywords internal
dist_to_diagonal <- function(X, ground = c("L2", "Linf")) {
  ground <- match.arg(ground)
  gap <- abs(X[, "death"] - X[, "birth"])
  if (ground == "L2") gap / sqrt(2) else gap / 2
}

#' Build the augmented (n+m) x (n+m) assignment cost matrix for two diagrams
#'
#' @param X,Y Matrices with columns \code{birth}, \code{death} (as produced
#'   by the internal \code{diagram_points()} helper).
#' @param ground Ground metric on the birth-death plane, \code{"L2"} or
#'   \code{"Linf"}.
#' @param power Exponent applied to every ground distance before it enters
#'   the matrix (\code{p} for a \eqn{p}-Wasserstein distance; use \code{1}
#'   for bottleneck, which works with raw distances and takes a max instead
#'   of a sum).
#' @return A square numeric matrix of size \code{nrow(X) + nrow(Y)}.
#'
#' @keywords internal
augmented_cost_matrix <- function(X, Y, ground = c("L2", "Linf"), power = 1) {
  ground <- match.arg(ground)
  n <- nrow(X); m <- nrow(Y)

  real_real <- if (n > 0 && m > 0) ground_dist(X, Y, ground) else matrix(0, n, m)
  dX <- if (n > 0) dist_to_diagonal(X, ground) else numeric(0)
  dY <- if (m > 0) dist_to_diagonal(Y, ground) else numeric(0)

  real_real <- real_real^power
  dX <- dX^power
  dY <- dY^power

  top <- cbind(real_real, matrix(dX, nrow = n, ncol = n))
  bottom <- cbind(matrix(dY, nrow = m, ncol = m, byrow = TRUE), matrix(0, nrow = m, ncol = n))
  rbind(top, bottom)
}

#' @keywords internal
wasserstein_assignment <- function(C, p) {
  assignment <- clue::solve_LSAP(C)
  total <- sum(C[cbind(seq_len(nrow(C)), assignment)])
  list(distance = total^(1 / p), assignment = assignment)
}

#' Wasserstein distance between two persistence diagrams
#'
#' The \eqn{p}-Wasserstein distance between the points of two persistence
#' diagrams in a given homological dimension, allowing points to be matched
#' to the diagonal (Cohen-Steiner, Edelsbrunner, Harer and Mileyko (2010),
#' "Lipschitz Functions Have \eqn{L_p}-Stable Persistence"). Essential
#' (death = \code{Inf}) classes are kept rather than dropped - see Details.
#' Computed exactly via the Hungarian algorithm (\code{clue::solve_LSAP}) on
#' the augmented assignment problem described in
#' \code{\link{augmented_cost_matrix}}.
#'
#' @param df1,df2 Persistence diagram data frames with columns \code{dim},
#'   \code{birth}, \code{death} (e.g. the output of
#'   \code{extract_persistence_pairs()} or \code{persistence_pairs()}).
#' @param dimension Homological dimension to compare.
#' @param p Wasserstein order (\eqn{p \ge 1}). Defaults to \code{2}.
#' @param ground Ground metric on the birth-death plane used to measure the
#'   distance between two points (or a point and the diagonal): \code{"L2"}
#'   (Euclidean, default) or \code{"Linf"} (Chebyshev).
#' @return A single non-negative number, the \eqn{p}-Wasserstein distance.
#'
#' @importFrom clue solve_LSAP
#' @export
#' @examples
#' df1 <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(1, 2))
#' df2 <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(1, 3))
#' wasserstein_distance(df1, df2, dimension = 0, p = 2)
wasserstein_distance <- function(df1, df2, dimension = 0, p = 2, ground = c("L2", "Linf")) {
  ground <- match.arg(ground)
  if (p < 1) stop("p must be >= 1.")

  X <- diagram_points(df1, dimension)
  Y <- diagram_points(df2, dimension)
  if (nrow(X) == 0 && nrow(Y) == 0) return(0)
  capped <- cap_essential(X, Y, dimension)
  X <- capped$X; Y <- capped$Y

  C <- augmented_cost_matrix(X, Y, ground = ground, power = p)
  wasserstein_assignment(C, p)$distance
}

#' @keywords internal
bipartite_perfect_matching_exists <- function(adj) {
  n <- nrow(adj)
  match_right <- integer(n) # 0 = unmatched

  augment <- function(u, visited) {
    for (v in which(adj[u, ])) {
      if (visited[v]) next
      visited[v] <- TRUE
      if (match_right[v] == 0L || augment(match_right[v], visited)) {
        match_right[v] <<- u
        return(TRUE)
      }
    }
    FALSE
  }

  for (u in seq_len(n)) {
    if (!augment(u, logical(n))) return(FALSE)
  }
  TRUE
}

#' @keywords internal
bipartite_maximum_matching <- function(adj) {
  # Same augmenting-path search as bipartite_perfect_matching_exists(), but
  # runs every row to completion instead of bailing out on the first row
  # that fails to match - used to recover the actual matching (not just
  # existence) once a feasible threshold has already been found.
  N <- nrow(adj)
  match_right <- integer(N) # 0 = unmatched; match_right[v] = matched row u

  augment <- function(u, visited) {
    for (v in which(adj[u, ])) {
      if (visited[v]) next
      visited[v] <- TRUE
      if (match_right[v] == 0L || augment(match_right[v], visited)) {
        match_right[v] <<- u
        return(TRUE)
      }
    }
    FALSE
  }

  for (u in seq_len(N)) augment(u, logical(N))
  match_right
}

#' @keywords internal
bottleneck_threshold <- function(C) {
  candidates <- sort(unique(c(0, as.vector(C))))
  lo <- 1L
  hi <- length(candidates)
  while (lo < hi) {
    mid <- (lo + hi) %/% 2L
    adj <- C <= candidates[mid] + 1e-9
    if (bipartite_perfect_matching_exists(adj)) hi <- mid else lo <- mid + 1L
  }
  candidates[lo]
}

#' Bottleneck distance between two persistence diagrams
#'
#' The bottleneck distance between the points of two persistence diagrams in
#' a given homological dimension, allowing points to be matched to the
#' diagonal: the infimum, over all matchings, of the largest single
#' point-to-point distance (Cohen-Steiner, Edelsbrunner and Harer (2007),
#' "Stability of Persistence Diagrams") - the \eqn{p \to \infty} limit of
#' \code{\link{wasserstein_distance}}. Essential (death = \code{Inf})
#' classes are kept rather than dropped - see
#' \code{\link{wasserstein_distance}}'s Details for why and how. Computed
#' exactly by binary search over the (finitely many) candidate distance
#' values for the smallest one admitting a perfect matching, checked with a
#' standard augmenting-path bipartite matcher, on the same augmented
#' assignment problem \code{\link{wasserstein_distance}} uses (see
#' \code{\link{augmented_cost_matrix}}).
#'
#' @param df1,df2 Persistence diagram data frames with columns \code{dim},
#'   \code{birth}, \code{death}.
#' @param dimension Homological dimension to compare.
#' @param ground Ground metric on the birth-death plane: \code{"L2"} or
#'   \code{"Linf"} (Chebyshev, the convention used by the reference above;
#'   default).
#' @return A single non-negative number, the bottleneck distance.
#'
#' @export
#' @examples
#' df1 <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(1, 2))
#' df2 <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(1, 3))
#' bottleneck_distance(df1, df2, dimension = 0)
bottleneck_distance <- function(df1, df2, dimension = 0, ground = c("Linf", "L2")) {
  ground <- match.arg(ground)

  X <- diagram_points(df1, dimension)
  Y <- diagram_points(df2, dimension)
  if (nrow(X) == 0 && nrow(Y) == 0) return(0)
  capped <- cap_essential(X, Y, dimension)
  X <- capped$X; Y <- capped$Y

  C <- augmented_cost_matrix(X, Y, ground = ground, power = 1)
  bottleneck_threshold(C)
}

#' Optimal point matching between two persistence diagrams
#'
#' Computes the same optimal matching that
#' \code{\link{wasserstein_distance}}/\code{\link{bottleneck_distance}}
#' reduce to internally, and returns it decoded into a readable table
#' instead of just the distance - what \code{\link{plot_matching}} draws.
#' Every point of both diagrams appears in exactly one row: matched to a
#' point of the other diagram (\code{type = "real-real"}), or matched to the
#' diagonal, i.e. effectively unmatched (\code{type = "x-diagonal"} for a
#' point of \code{df1}, \code{"y-diagonal"} for a point of \code{df2}).
#'
#' @param df1,df2 Persistence diagram data frames with columns \code{dim},
#'   \code{birth}, \code{death}.
#' @param dimension Homological dimension to compare.
#' @param distance Which distance's optimal matching to compute,
#'   \code{"wasserstein"} (default) or \code{"bottleneck"}.
#' @param p Wasserstein order, used only when \code{distance = "wasserstein"}.
#' @param ground Ground metric on the birth-death plane, \code{"L2"} or
#'   \code{"Linf"}. Defaults to each distance's own default (\code{"L2"} for
#'   Wasserstein, \code{"Linf"} for bottleneck) when left \code{NULL}.
#' @return A list with \code{distance} (the matching's cost, equal to what
#'   \code{\link{wasserstein_distance}}/\code{\link{bottleneck_distance}}
#'   would return), \code{matches} (a data frame with columns
#'   \code{x_birth, x_death, x_essential, y_birth, y_death, y_essential,
#'   type}; a \code{NA} pair on one side means that row's point matched the
#'   diagonal), and \code{X}, \code{Y}, \code{essential_X}, \code{essential_Y}
#'   (the two diagrams' points after the same essential-pair capping
#'   \code{\link{wasserstein_distance}} uses - see its Details - and which of
#'   them were essential before capping).
#'
#' @export
#' @examples
#' df1 <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(1, 2))
#' df2 <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(1, 3))
#' diagram_matching(df1, df2, dimension = 0, distance = "wasserstein", p = 2)
diagram_matching <- function(df1, df2, dimension, distance = c("wasserstein", "bottleneck"),
                              p = 2, ground = NULL) {
  distance <- match.arg(distance)
  if (is.null(ground)) {
    ground <- if (distance == "bottleneck") "Linf" else "L2"
  } else {
    ground <- match.arg(ground, c("L2", "Linf"))
  }
  if (distance == "wasserstein" && p < 1) stop("p must be >= 1.")

  X <- diagram_points(df1, dimension)
  Y <- diagram_points(df2, dimension)
  n <- nrow(X); m <- nrow(Y)
  if (n == 0 && m == 0) {
    stop(sprintf("No points in dimension %d in either diagram - nothing to match.", dimension))
  }

  essential_X <- if (n > 0) !is.finite(X[, "death"]) else logical(0)
  essential_Y <- if (m > 0) !is.finite(Y[, "death"]) else logical(0)

  capped <- cap_essential(X, Y, dimension)
  X <- capped$X; Y <- capped$Y

  power <- if (distance == "wasserstein") p else 1
  C <- augmented_cost_matrix(X, Y, ground = ground, power = power)

  if (distance == "wasserstein") {
    solved <- wasserstein_assignment(C, p)
    dist_val <- solved$distance
    row_to_col <- solved$assignment # row_to_col[i] = column matched to row i
  } else {
    dist_val <- bottleneck_threshold(C)
    adj <- C <= dist_val + 1e-9
    match_right <- bipartite_maximum_matching(adj) # match_right[v] = row matched to column v
    row_to_col <- integer(length(match_right))
    row_to_col[match_right] <- seq_along(match_right)
  }

  # Decode: rows 1..n are df1's real points, rows n+1..n+m are diagonal
  # slots that let df2's real points match the diagonal; columns 1..m are
  # df2's real points, columns m+1..m+n are diagonal slots for df1's real
  # points (see augmented_cost_matrix()).
  matches <- vector("list", n + m)
  k <- 0L
  for (i in seq_len(n)) {
    j <- row_to_col[i]
    k <- k + 1L
    if (j <= m) {
      matches[[k]] <- data.frame(
        x_birth = X[i, "birth"], x_death = X[i, "death"], x_essential = essential_X[i],
        y_birth = Y[j, "birth"], y_death = Y[j, "death"], y_essential = essential_Y[j],
        type = "real-real"
      )
    } else {
      matches[[k]] <- data.frame(
        x_birth = X[i, "birth"], x_death = X[i, "death"], x_essential = essential_X[i],
        y_birth = NA_real_, y_death = NA_real_, y_essential = NA,
        type = "x-diagonal"
      )
    }
  }
  if (m > 0) {
    for (v in seq_len(m)) {
      row_idx <- n + v
      j <- row_to_col[row_idx]
      if (j <= m) {
        # a diagonal-slot row matched to column j means df2's point j went
        # to the diagonal; the row's own identity (v) is not meaningful
        k <- k + 1L
        matches[[k]] <- data.frame(
          x_birth = NA_real_, x_death = NA_real_, x_essential = NA,
          y_birth = Y[j, "birth"], y_death = Y[j, "death"], y_essential = essential_Y[j],
          type = "y-diagonal"
        )
      }
      # else: this diagonal slot matched a df1 diagonal slot (cost 0,
      # neither side a real point) - nothing to report
    }
  }

  list(
    distance = dist_val,
    matches = do.call(rbind, matches[seq_len(k)]),
    X = X, Y = Y, essential_X = essential_X, essential_Y = essential_Y
  )
}
