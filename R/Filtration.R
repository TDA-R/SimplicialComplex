#' Vietoris-Rips Filtration: Get the boundary matrix and its reduction information in matrix form
#'
#' @param points Data point input.
#' @param eps_max Maximum scale (epsilon).
#' @return A list of simplices with their information.
#'
#' @export
#' @examples
#' points <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1), ncol = 2)
#' filtration <- build_vr_filtration(points, eps_max=1.2)
build_vr_filtration <- function(points, eps_max) {
  vr <- VietorisRipsComplex(points, epsilon = eps_max)
  kmax <- max(sapply(vr$simplices, length)) - 1 # Dim of each simplices
  face <- lapply(0:kmax, function(k) faces(vr$simplices, k))
  all_faces <- unlist(face,recursive = FALSE)

  # Give information to each simplex
  simplex_info <- lapply(all_faces, function(s) {
    list(simplex=s, t=vr_scale_of_simplex(points, s))
  })
  ord <- order(sapply(simplex_info, function(x) x[["t"]]),
               sapply(simplex_info, function(x) length(x$simplex)),
               sapply(simplex_info, function(x) paste(x$simplex, collapse="-")))

  return(simplex_info[ord])
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



