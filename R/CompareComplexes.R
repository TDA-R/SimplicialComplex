#' Compare simplicial complex constructions on one point cloud
#'
#' Builds the same point cloud's filtration under several
#' \code{\link{build_filtration}} methods and summarizes, side by side, how
#' expensive each one was to build and what persistent homology it found -
#' useful for seeing directly how Cech's dimension blow-up, Alpha's
#' Delaunay-bounded dimension, VR's cheap-but-approximate cliques, and
#' Witness's landmark subsampling trade off against each other on the same
#' data (see Otter et al. (2017), Section 5.2, for the underlying
#' trade-offs).
#'
#' @param points A numeric matrix or data.frame, one point per row.
#' @param methods Character vector of methods to compare, any subset of
#'   \code{"VR"}, \code{"Delaunay"}, \code{"Alpha"}, \code{"Cech"},
#'   \code{"Witness"}. Defaults to all five.
#' @param eps_max Maximum scale, passed to \code{\link{build_filtration}}.
#'   Required whenever \code{methods} includes \code{"VR"}, \code{"Cech"} or
#'   \code{"Witness"}; optional for \code{"Alpha"} (defaults to \code{Inf});
#'   ignored for \code{"Delaunay"}.
#' @param landmarks Passed to \code{\link{build_filtration}} for
#'   \code{"Witness"}. If \code{NULL} and \code{"Witness"} is included,
#'   defaults to \code{min(30, nrow(points))} landmarks via
#'   \code{\link{generate_landmarks}}, with a message.
#' @param nu Passed to \code{\link{build_filtration}} for \code{"Witness"}.
#' @return A list with:
#' \describe{
#'   \item{summary}{A data frame, one row per method, with the number of
#'     vertices and simplices, the highest simplex dimension reached, build
#'     and persistence-computation time in seconds, and the number of
#'     finite/essential persistence pairs found.}
#'   \item{diagrams}{A named list of persistence diagram data frames (one per
#'     method, from \code{\link{persistence_pairs}}), for further comparison
#'     (e.g. with \code{\link{wasserstein_distance}}/
#'     \code{\link{bottleneck_distance}}, or \code{\link{plot_persistence}}).}
#'   \item{filtrations}{A named list of the raw filtration lists.}
#' }
#'
#' @export
#' @examples
#' points <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1, 0.5, 0.5), ncol = 2, byrow = TRUE)
#' cmp <- compare_complexes(points, methods = c("VR", "Alpha", "Delaunay"), eps_max = 1.2)
#' cmp$summary
compare_complexes <- function(points,
                               methods = c("VR", "Delaunay", "Alpha", "Cech", "Witness"),
                               eps_max = NULL,
                               landmarks = NULL,
                               nu = 1) {
  points <- as.matrix(points)
  methods <- match.arg(methods,
                        choices = c("VR", "Delaunay", "Alpha", "Cech", "Witness"),
                        several.ok = TRUE)

  if (any(methods %in% c("VR", "Cech", "Witness")) && is.null(eps_max)) {
    stop("eps_max is required for methods 'VR', 'Cech' and 'Witness'.")
  }
  if ("Witness" %in% methods && is.null(landmarks)) {
    landmarks <- min(30L, nrow(points))
    message(sprintf(
      "No landmarks supplied for method = 'Witness'; using %d landmarks via generate_landmarks().",
      landmarks
    ))
  }

  filtrations <- list()
  diagrams <- list()
  rows <- vector("list", length(methods))

  for (i in seq_along(methods)) {
    m <- methods[i]
    build_time <- system.time(
      filt <- if (m == "Delaunay") {
        build_filtration(points, method = m)
      } else if (m == "Witness") {
        build_filtration(points, method = m, eps_max = eps_max, landmarks = landmarks, nu = nu)
      } else {
        build_filtration(points, method = m, eps_max = eps_max)
      }
    )[["elapsed"]]

    pers_time <- system.time(diag <- persistence_pairs(filt))[["elapsed"]]

    dims <- lengths(lapply(filt, `[[`, "simplex")) - 1L

    filtrations[[m]] <- filt
    diagrams[[m]] <- diag
    rows[[i]] <- data.frame(
      method = m,
      n_vertices = sum(dims == 0),
      n_simplices = length(filt),
      max_dimension = max(dims),
      build_time_sec = round(build_time, 4),
      persistence_time_sec = round(pers_time, 4),
      n_finite_pairs = sum(is.finite(diag$death)),
      n_essential_pairs = sum(is.infinite(diag$death)),
      stringsAsFactors = FALSE
    )
  }

  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL

  list(summary = summary, diagrams = diagrams, filtrations = filtrations)
}
