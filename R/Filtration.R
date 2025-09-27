library(igraph)

source('./R/VRComplex.R')
source('./R/Faces.R')

vr_scale_of_simplex <- function(
    points, simplex_idx
    ) {
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
  dmax
}

build_vr_filtration <- function(
    points, eps_max
    ) {
  vr <- VietorisRipsComplex(points, epsilon = eps_max)
  kmax <- max(sapply(vr$simplices, length)) - 1 # Dim of each simplices
  face <- lapply(0:kmax, function(k) faces(vr$simplices, k))
  all_faces <- unlist(face,recursive = FALSE)

  # Give information to each simplex
  F <- lapply(all_faces, function(s) {
    list(simplex=s, t=vr_scale_of_simplex(points, s))
  })
  ord <- order(sapply(F, function(x) x[["t"]]),
               sapply(F, function(x) length(x$simplex)),
               sapply(F, function(x) paste(x$simplex, collapse="-")))
  F <- F[ord]
}

# 15 points
points <- matrix(c(
  0,0, 1,0, 0,1, 1,1,
  0.5,0.5, 0.5,1.5, 1.5,0.5, 1.5,1.5,
  3,3, 4,3, 3,4, 4,4
), ncol=2, byrow=TRUE)

F <- build_vr_filtration(points, eps_max=1.2)



