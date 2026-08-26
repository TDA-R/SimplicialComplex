#' Wasserstein/bottleneck distance between two point clouds, via one chosen complex
#'
#' Picks a single \code{\link{build_filtration}} method, builds the
#' persistence diagram of each point cloud with it, and compares the two
#' diagrams with either \code{\link{wasserstein_distance}} or
#' \code{\link{bottleneck_distance}}. Useful for asking "how different are
#' these two datasets topologically", holding the complex construction fixed
#' so the comparison is apples-to-apples.
#'
#' @param points1,points2 Numeric matrices or data.frames, one point per row.
#' @param method One of \code{"VR"}, \code{"Delaunay"}, \code{"Alpha"},
#'   \code{"Cech"}, \code{"Witness"} - passed to \code{\link{build_filtration}}.
#' @param distance One of \code{"wasserstein"} or \code{"bottleneck"}.
#' @param dimension Homological dimension to compare.
#' @param p Wasserstein order, used only when \code{distance = "wasserstein"}.
#' @param ground Ground metric passed to the chosen distance function
#'   (\code{"L2"} or \code{"Linf"}). Defaults to each distance's own default
#'   (\code{"L2"} for Wasserstein, \code{"Linf"} for bottleneck) if left
#'   \code{NULL}.
#' @param eps_max Maximum scale, passed to \code{\link{build_filtration}}.
#'   Required for \code{"VR"}, \code{"Cech"} and \code{"Witness"}; optional
#'   for \code{"Alpha"} (defaults to \code{Inf}); ignored for
#'   \code{"Delaunay"}.
#' @param landmarks Passed to \code{\link{build_filtration}} for
#'   \code{"Witness"}; if \code{NULL}, defaults independently for each point
#'   cloud to \code{min(30, nrow(points))} landmarks, with a message.
#' @param nu Passed to \code{\link{build_filtration}} for \code{"Witness"}.
#' @return A list with:
#' \describe{
#'   \item{distance}{The computed distance (a single number).}
#'   \item{method, distance_type, dimension}{Echoed back for reference.}
#'   \item{diagram1, diagram2}{The two persistence diagrams (data frames)
#'     the distance was computed from.}
#' }
#'
#' @export
#' @examples
#' set.seed(1)
#' cloud_a <- matrix(rnorm(24), ncol = 2)
#' cloud_b <- matrix(rnorm(24), ncol = 2) + 0.1
#' complex_distance(cloud_a, cloud_b, method = "VR", distance = "bottleneck",
#'                   dimension = 0, eps_max = 0.6)
complex_distance <- function(points1, points2,
                              method = c("VR", "Delaunay", "Alpha", "Cech", "Witness"),
                              distance = c("wasserstein", "bottleneck"),
                              dimension = 0,
                              p = 2,
                              ground = NULL,
                              eps_max = NULL,
                              landmarks = NULL,
                              nu = 1) {
  method <- match.arg(method)
  distance <- match.arg(distance)
  points1 <- as.matrix(points1)
  points2 <- as.matrix(points2)

  if (method %in% c("VR", "Cech", "Witness") && is.null(eps_max)) {
    stop(sprintf("eps_max is required for method = '%s'.", method))
  }

  build_one <- function(points, lm) {
    if (method == "Delaunay") {
      build_filtration(points, method = method)
    } else if (method == "Witness") {
      build_filtration(points, method = method, eps_max = eps_max, landmarks = lm, nu = nu)
    } else {
      build_filtration(points, method = method, eps_max = eps_max)
    }
  }

  lm1 <- landmarks
  lm2 <- landmarks
  if (method == "Witness" && is.null(landmarks)) {
    lm1 <- min(30L, nrow(points1))
    lm2 <- min(30L, nrow(points2))
    message(sprintf(
      "No landmarks supplied for method = 'Witness'; using %d and %d landmarks (points1, points2 respectively).",
      lm1, lm2
    ))
  }

  diagram1 <- persistence_pairs(build_one(points1, lm1))
  diagram2 <- persistence_pairs(build_one(points2, lm2))

  if (is.null(ground)) ground <- if (distance == "bottleneck") "Linf" else "L2"

  value <- if (distance == "wasserstein") {
    wasserstein_distance(diagram1, diagram2, dimension = dimension, p = p, ground = ground)
  } else {
    bottleneck_distance(diagram1, diagram2, dimension = dimension, ground = ground)
  }

  list(distance = value, method = method, distance_type = distance, dimension = dimension,
       diagram1 = diagram1, diagram2 = diagram2)
}
