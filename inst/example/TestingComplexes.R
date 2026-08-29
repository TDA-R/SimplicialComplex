library(Matrix)
library(gtools)
library(igraph)
library(ggplot2)
library(geometry)
library(clue)

source("./R/Faces.R")
source("./R/Betti.R")
source("./R/Boundary.R")
source("./R/EulerCharacteristic.R")
source("./R/ComplexUtils.R")
source("./R/AlphaComplex.R")
source("./R/CechComplex.R")
source("./R/FloodComplex.R") # generate_landmarks(), as_filtration()
source("./R/WitnessComplex.R")
source("./R/FloodComplex.R")
source("./R/Filtration.R") # build_filtration()
source("./R/Persistence.R")
source("./R/PlotPD.R")

set.seed(42)
points <- matrix(runif(20), ncol = 2)

# Vietoris-Rips (existing)
vr_filtration <- build_filtration(points, method = "VR", eps_max = 1.2, max_dimension=2)
cat(sprintf("VR: %d simplices\n", length(vr_filtration)))

# Delaunay
delaunay_filtration <- build_filtration(points, method = "Delaunay", max_dimension=2)
cat(sprintf("Delaunay: %d simplices\n", length(delaunay_filtration)))

# Alpha
alpha_filtration <- build_filtration(points, method = "Alpha", eps_max = 1.2, max_dimension=2)
cat(sprintf("Alpha: %d simplices\n", length(alpha_filtration)))

# Cech
cech_filtration <- build_filtration(points, method = "Cech", eps_max = 1.2, max_dimension=2)
cech_max_dim <- max(lengths(lapply(cech_filtration, `[[`, "simplex")) - 1L)
cat(sprintf("Cech: %d simplices (max dimension %d)\n", length(cech_filtration), cech_max_dim))

# Witness
theta <- seq(0, 2 * pi, length.out = 40)[-40]
ring <- cbind(cos(theta), sin(theta)) + matrix(rnorm(2 * length(theta), sd = 0.03), ncol = 2)
n_landmarks <- 10
witness_filtration <- build_filtration(ring, method = "Witness", eps_max = 0.6, max_dimension=2, landmarks = n_landmarks)
cat(sprintf("Witness: %d simplices on %d landmarks\n", length(witness_filtration), n_landmarks))

# Flood
flood_filtration <- build_flood_filtration(points, landmarks = 20, points_per_edge = 15, backend = "cpu")
cat(sprintf("Flood: %d simplices\n", length(flood_filtration)))


show_top_bars <- function(pairs, label, top_n = 3) {
  pairs$pers <- pairs$death - pairs$birth
  cat(sprintf("\n--- %s ---\n", label))
  for (d in sort(unique(pairs$dim))) {
    pd <- pairs[pairs$dim == d, ]
    pd <- pd[order(-pd$pers), ]
    cat(sprintf("H%d (top %d of %d):\n", d, min(top_n, nrow(pd)), nrow(pd)))
    print(head(pd[, c("birth", "death", "pers")], top_n), row.names = FALSE)
  }
}

show_top_bars(persistence_pairs(vr_filtration), "VR")
show_top_bars(persistence_pairs(delaunay_filtration), "Delaunay")
show_top_bars(persistence_pairs(alpha_filtration), "Alpha")
show_top_bars(persistence_pairs(cech_filtration), "Cech")
show_top_bars(persistence_pairs(witness_filtration), "Witness")

plot_persistence(persistence_pairs(vr_filtration))

# Cubical
source("R/CubicalComplex.R")
image <- matrix(c(
  9, 9, 9, 9, 9,
  9, 1, 1, 1, 9,
  9, 1, 5, 1, 9,
  9, 1, 1, 1, 9,
  9, 9, 9, 9, 9
), nrow = 5, byrow = TRUE)

filtration <- build_cubical_filtration(image)
pairs <- persistence_pairs(filtration)
pairs[pairs$dim == 1, ]
plot_persistence(pairs)
show_top_bars(pairs, "Cubical")
