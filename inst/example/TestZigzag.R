library(Matrix)
library(gtools)

source("R/Faces.R")
source("R/Boundary.R")
source("R/Betti.R")
source("R/Persistence.R")
source("R/zigzag.R")

# Normal example
K0 <- list(c(1, 2), c(2, 3), c(3, 4), c(4, 1), c(1, 3)) # square + diagonal: H1 = 2
K1 <- list(c(1, 2), c(2, 3), c(3, 4), c(4, 1))          # delete diagonal:   H1 = 1
K2 <- list(c(1, 2), c(2, 3), c(3, 4), c(4, 1), c(1, 3)) # re-insert diagonal: H1 = 2
K3 <- list(c(1, 2), c(2, 3), c(3, 4), c(1, 3))          # delete edge 4-1:    H1 = 1
K4 <- list(c(1, 2), c(2, 3), c(1, 3), 4)                # delete edge 3-4 (bridge); H0 splits

Ks <- list(K0, K1, K2, K3, K4)
bars <- zigzag_persistence(Ks)
print(bars)

plot_persistence(bars)


# Purely-growing example
K0 <- list(1)
K1 <- list(1, 2)
K2 <- list(1, 2, 3)
K3 <- list(1, 2, 3, c(1, 2))
K4 <- list(1, 2, 3, c(1, 2), c(1, 3))
K5 <- list(1, 2, 3, c(1, 2), c(1, 3), c(2, 3))
K6 <- list(1, 2, 3, c(1, 2), c(1, 3), c(2, 3), c(1, 2, 3))
n <- 6

zz <- zigzag_persistence(list(K0, K1, K2, K3, K4, K5, K6))
print(zz)

plot_persistence(zz)
