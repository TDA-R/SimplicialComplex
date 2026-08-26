#' Build a filtration from a point cloud, for any of five complex types
#'
#' @param points A numeric matrix or data.frame, one point per row.
#' @param method One of \code{"VR"}, \code{"Delaunay"}, \code{"Alpha"},
#'   \code{"Cech"}, \code{"Witness"}.
#' @param eps_max Maximum scale. Required for \code{"VR"}, \code{"Cech"} and
#'   \code{"Witness"}; for \code{"Alpha"} it caps the alpha value (defaults
#'   to \code{Inf}, i.e. no cap); ignored for \code{"Delaunay"} (which has no
#'   scale parameter - see \code{\link{DelaunayComplex}}).
#' @param landmarks For \code{method = "Witness"} only: either an integer
#'   (number of landmarks to select via \code{\link{generate_landmarks}}) or
#'   an integer vector of landmark row indices into \code{points}. Required
#'   for \code{"Witness"}.
#' @param nu For \code{method = "Witness"} only: the landmark-distance
#'   parameter passed to \code{\link{WitnessComplex}}. Defaults to \code{1}.
#' @param max_dimension Optional maximum homology dimension you intend to
#'   compute persistence for.
#' @return A filtration list, as described above. When \code{max_dimension}
#'   is set, the list carries it as a \code{"max_dimension"} attribute (see
#'   above); otherwise the attribute is absent, exactly as before this
#'   parameter existed.
#'
#' @export
#' @examples
#' points <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1), ncol = 2)
#' vr_filtration <- build_filtration(points, method = "VR", eps_max = 1.2)
#' alpha_filtration <- build_filtration(points, method = "Alpha")
#' delaunay_filtration <- build_filtration(points, method = "Delaunay")
#' cech_filtration <- build_filtration(points, method = "Cech", eps_max = 0.8)
#' # cap at H1: persistence_pairs()/etc. need no max_dimension of their own
#' vr_capped <- build_filtration(points, method = "VR", eps_max = 1.2,
#'                                max_dimension = 1)
#' pairs <- persistence_pairs(vr_capped)
build_filtration <- function(points, method = c("VR", "Delaunay", "Alpha", "Cech", "Witness"),
                              eps_max = NULL, landmarks = NULL, nu = 1, max_dimension = NULL) {
  method <- match.arg(method)
  points <- as.matrix(points)

  if (method %in% c("VR", "Cech", "Witness") && is.null(eps_max)) {
    stop(sprintf("method = '%s' requires eps_max.", method))
  }
  if (method == "Witness" && is.null(landmarks)) {
    stop("method = 'Witness' requires landmarks (a count or a vector of row indices).")
  }

  # The "+1 trick": one extra dimension kept internally so dimension-
  # max_dimension classes are still correctly killed - see
  # restrict_filtration()'s Details. NULL (uncapped) stays NULL.
  build_kmax <- if (is.null(max_dimension)) NULL else max_dimension + 1

  result <- switch(
    method,
    VR = {
      vr <- VietorisRipsComplex(points, epsilon = eps_max)
      simplices_to_filtration(vr$simplices, function(s) vr_scale_of_simplex(points, s),
                               max_dimension = build_kmax)
    },
    Cech = {
      cc <- CechComplex(points, epsilon = eps_max)
      simplices_to_filtration(cc$simplices, function(s) cech_scale_of_simplex(points, s),
                               max_dimension = build_kmax)
    },
    Witness = {
      wc <- WitnessComplex(points, landmarks = landmarks, epsilon = eps_max, nu = nu)
      simplices_to_filtration(wc$simplices, function(s) witness_scale_of_simplex(wc$edge_birth, s),
                               max_dimension = build_kmax)
    },
    Alpha = {
      epsilon <- if (is.null(eps_max)) Inf else eps_max
      alpha_filtration <- as_filtration(AlphaComplex(points, epsilon = epsilon))
      # AlphaComplex() has no build-time dimension cap, so restrict after
      # the fact instead - still correct (just not cheaper to build).
      if (is.null(build_kmax)) alpha_filtration else restrict_filtration(alpha_filtration, build_kmax)
    },
    Delaunay = {
      delaunay_filtration <- as_filtration(DelaunayComplex(points))
      if (is.null(build_kmax)) delaunay_filtration else restrict_filtration(delaunay_filtration, build_kmax)
    }
  )

  # restrict_filtration() filters a plain list (via Filter()/subsetting),
  # which drops attributes - so max_dimension must be tagged here, after
  # the switch, not earlier. NULL is a safe no-op when no cap was requested.
  attr(result, "max_dimension") <- max_dimension
  result
}

#' @keywords internal
vr_scale_of_simplex <- function(points, simplex_idx) {
  # Get the maximum distance between points in the simplex
  if (length(simplex_idx) <= 1) return(0)
  # Get specific rows
  srow <- points[simplex_idx, , drop=FALSE]
  dmax <- 0
  # Loop will give all the conditions
  # If nrow = 4, 1 -> 2,3,4. 2 -> 3,4. 3 -> 4
  for (i in 1:(nrow(srow)-1)) {
    for (j in (i+1):nrow(srow)) {
      d <- sqrt(sum((srow[i,]-srow[j,])^2))
      if (d > dmax) dmax <- d
    }
  }
  return(dmax)
}
