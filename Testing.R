library(Matrix)
library(gtools)
library(igraph)

source("./R/Faces.R")
source("./R/Betti.R")
source("./R/Boundary.R")
source("./R/EulerCharacteristic.R")
source("./R/AbstractSimplicialComplex.R")

simplices <- list(c(1, 2), c(3, 4), c(2, 1, 3), c(4, 2))

faces(simplices, target_dim=0)
faces(simplices, target_dim=1)
faces(simplices, target_dim=2)
faces(simplices, target_dim=3)

boundary(simplices, 0)
boundary(simplices, 1)
boundary(simplices, 2)

betti_number(simplices, 0, tol=0.1)
betti_number(simplices, 1, tol=0.1)
betti_number(simplices, 2, tol=0.1)
betti_number(simplices, 3, tol=0.1)

euler_characteristic(simplices, tol=0.1)

abstract_simplicial_complex(simplices, 2)

source("R/VRComplex.R")

points <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1), ncol = 2)
epsilon <- 1.5

vr_complex <- VietorisRipsComplex(points, epsilon)

# E(vr_complex$network)
# vr_complex$simplices

# Topological Structure
faces(vr_complex$simplices, target_dim=2)
boundary(vr_complex$simplices, 2)
betti_number(vr_complex$simplices, 1, tol=0.1)

plot(
  vr_complex$network,
  layout = points,
  vertex.label = 1:nrow(points),
  vertex.size = 12,
  edge.arrow.mode = 0,
  asp = 1
)
