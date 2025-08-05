source("R/faces.R")
source("R/Betti.R")
source("R/Boundary.R")
source("R/EulerCharacteristic")
source("R/AbstractSimplicialComplex.R")

simplices <- list(c(2, 1, 3), c(4, 2), c(5), c(2, 3, 5, 4))

faces(simplices, target_dim=0)
faces(simplices, target_dim=1)
faces(simplices, target_dim=2)
faces(simplices, target_dim=3)

boundary(simplices, 0)
boundary(simplices, 1)
boundary(simplices, 2)

betti_number(simplices, 0, eps=0.1)
betti_number(simplices, 1, eps=0.1)
betti_number(simplices, 2, eps=0.1)
betti_number(simplices, 3, eps=0.1)

euler_characteristic(simplices, eps=0.1)

abstract_simplicial_complex(simplices, 2)
