#' Build a cubical (grid) filtration from an image
#'
#' @param image A numeric matrix (grid of pixel/voxel values, e.g. a
#'   grayscale image with values in \code{[0, 255]}).
#' @param superlevel If \code{TRUE}, filters by decreasing value instead
#'   (equivalent to negating \code{image} first) - useful for tracking
#'   bright structures shrinking rather than dark structures growing.
#' @return A filtration list: one \code{list(simplex = <integer vector of
#'   pixel ids>, t = <numeric>)} per simplex, sorted by (t, dimension,
#'   lexicographic order) - the same format as
#'   \code{\link{build_filtration}}.
#'
#' @export
#' @examples
#' # a ring (value 1) around a hole (value 5) around a background (value 9):
#' # the hole is born once the ring closes and dies once its center fills in
#' image <- matrix(c(
#'   9, 9, 9, 9, 9,
#'   9, 1, 1, 1, 9,
#'   9, 1, 5, 1, 9,
#'   9, 1, 1, 1, 9,
#'   9, 9, 9, 9, 9
#' ), nrow = 5, byrow = TRUE)
#' filtration <- build_cubical_filtration(image)
#' pairs <- persistence_pairs(filtration)
#' pairs[pairs$dim == 1, ] # one H1 bar: birth = 1 (ring closes), death = 5
#' plot_persistence(pairs)
build_cubical_filtration <- function(image, superlevel = FALSE) {
  image <- as.matrix(image)
  nr <- nrow(image)
  nc <- ncol(image)

  value <- if (superlevel) -image else image

  vid <- function(i, j) (i - 1L) * nc + j # 1-based linear pixel id, row-major
  coord <- function(v) {
    v0 <- v - 1L
    c(v0 %/% nc + 1L, v0 %% nc + 1L) # c(row, col)
  }
  cell_value <- function(verts) {
    rc <- vapply(verts, coord, numeric(2))
    max(value[cbind(rc[1, ], rc[2, ])])
  }

  simplex_info <- vector("list", 5L * nr * nc)
  idx <- 0L
  for (i in seq_len(nr)) {
    for (j in seq_len(nc)) {
      idx <- idx + 1L
      simplex_info[[idx]] <- list(simplex = vid(i, j), t = value[i, j])
    }
  }

  seen <- new.env(hash = TRUE, parent = emptyenv())

  add_simplex <- function(verts) {
    key <- paste(sort(verts), collapse = "-")
    if (exists(key, envir = seen, inherits = FALSE)) return(invisible(NULL))
    assign(key, TRUE, envir = seen)
    idx <<- idx + 1L
    simplex_info[[idx]] <<- list(simplex = sort(verts), t = cell_value(verts))
  }

  for (i in seq_len(nr - 1L)) {
    for (j in seq_len(nc - 1L)) {
      v00 <- vid(i, j);      v10 <- vid(i + 1L, j)
      v01 <- vid(i, j + 1L); v11 <- vid(i + 1L, j + 1L)

      add_simplex(c(v00, v10)) # left edge
      add_simplex(c(v00, v01)) # top edge
      add_simplex(c(v10, v11)) # right edge
      add_simplex(c(v01, v11)) # bottom edge
      add_simplex(c(v00, v11)) # diagonal

      add_simplex(c(v00, v10, v11))
      add_simplex(c(v00, v01, v11))
    }
  }

  simplex_info <- simplex_info[seq_len(idx)]

  ord <- order(
    vapply(simplex_info, `[[`, 0, "t"),
    lengths(lapply(simplex_info, `[[`, "simplex")),
    vapply(simplex_info, function(x) paste(x$simplex, collapse = "-"), "")
  )
  simplex_info[ord]
}
