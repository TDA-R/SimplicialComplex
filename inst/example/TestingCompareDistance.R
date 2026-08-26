library(Matrix)
library(gtools)
library(igraph)
library(ggplot2)
library(geometry)
library(clue)

source("./R/Faces.R")
source("./R/Boundary.R")
source("./R/Betti.R")
source("./R/ComplexUtils.R")
source("./R/VRComplex.R")
source("./R/AlphaComplex.R")
source("./R/CechComplex.R")
source("./R/FloodComplex.R") # generate_landmarks(), as_filtration()
source("./R/WitnessComplex.R")
source("./R/Filtration.R") # build_filtration()
source("./R/Persistence.R")
source("./R/Distance.R") # wasserstein_distance(), bottleneck_distance()
source("./R/CompareComplexes.R")
source("./R/ComplexDistance.R")
source("./R/PlotPD.R")
source("./R/PlotMatching.R") # diagram_matching(), plot_matching()

set.seed(11)
theta <- seq(0, 2 * pi, length.out = 30)[-30]
ring <- cbind(cos(theta), sin(theta)) + matrix(rnorm(2 * length(theta), sd = 0.03), ncol = 2)
plot(ring)

cmp <- compare_complexes(ring, methods = c("VR", "Delaunay", "Alpha", "Cech"), eps_max = 0.5)
print(cmp$summary)

cmp_w <- compare_complexes(ring, methods = "Witness", eps_max = 0.5, landmarks = 10)
print(cmp_w$summary)

# most persistent H1 bar found by each method - essential (death = Inf)
# classes are kept, not filtered, same convention as TestingComplexes.R
for (m in names(cmp$diagrams)) {
  pd <- cmp$diagrams[[m]]
  pd <- pd[pd$dim == 1, ]
  pd$pers <- pd$death - pd$birth
  best <- pd[which.max(pd$pers), ]
  cat(sprintf("%-9s longest H1 bar: birth=%.3f death=%s\n", m,
              best$birth, ifelse(is.infinite(best$death), "Inf", sprintf("%.3f", best$death))))
}


ring_noisy <- cbind(cos(theta), sin(theta)) + matrix(rnorm(2 * length(theta), sd = 0.15), ncol = 2)
plot(ring_noisy)

res_w <- complex_distance(ring, ring_noisy, method = "VR", distance = "wasserstein", dimension = 1, p = 2, eps_max = 0.6)
cat(sprintf("\nVR Wasserstein-2 (H1), clean vs. noisy ring: %.4f\n", res_w$distance))

res_b <- complex_distance(ring, ring_noisy, method = "Alpha", distance = "bottleneck", dimension = 1)
cat(sprintf("Alpha bottleneck (H1), clean vs. noisy ring: %.4f\n", res_b$distance))

# self-distance sanity check: identical point clouds are always distance 0
res_self <- complex_distance(ring, ring, method = "Delaunay", distance = "wasserstein", dimension = 1, p = 2)
cat(sprintf("Delaunay self-distance (should be 0): %.4f\n", res_self$distance))

plot_matching(res_w$diagram1, res_w$diagram2, dimension = 1,
              distance = "wasserstein", p = 2,
              labels = c("Clean ring", "Noisy ring"))
