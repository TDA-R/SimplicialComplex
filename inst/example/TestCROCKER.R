library(Matrix)
library(gtools)
library(igraph)
library(ggplot2)
library(parallel)

source("R/Faces.R")
source("R/Boundary.R")
source("R/Betti.R")
source("R/ComplexUtils.R")
source("R/Filtration.R")
source("R/VRComplex.R")
source("R/Persistence.R")
source("R/Crocker.R")

# Scenario 1: a cluster of points dispersing over time
set.seed(42)
n_points <- 20
n_frames <- 15
base <- matrix(rnorm(n_points * 2, sd = 0.15), ncol = 2)
frames_dispersing <- lapply(seq_len(n_frames), function(t) {
  drift <- (t - 1) * 0.05
  jitter <- matrix(rnorm(n_points * 2, sd = 0.03), ncol = 2)
  outward <- matrix(rep(c(1, -1), length.out = n_points * 2), ncol = 2) * drift * 0.3
  base + jitter * (1 + drift) + outward
})

cr0 <- crocker(frames_dispersing, dim = 0, n_eps = 30, eps_max = 0.4)
plot_crocker(cr0) + labs(subtitle = "b0: dispersing point cluster")

# Scenario 2: a loop that closes then reopens (exercises dim = 1)
ring <- function(radius, n = 20) {
  theta <- seq(0, 2 * pi, length.out = n + 1)[1:n]
  cbind(radius * cos(theta), radius * sin(theta))
}
radii <- c(seq(1, 0.15, length.out = 6), seq(0.15, 1, length.out = 6))
frames_ring <- lapply(radii, ring)

cr1 <- crocker(frames_ring, dim = 1, n_eps = 30)
plot_crocker(cr1) + labs(subtitle = "b1: ring closing and reopening")

# Parallel vs sequential
n_cores_available <- parallel::detectCores()
cat("Cores detected:", n_cores_available, "\n")

t_seq <- system.time(cr_seq <- crocker(frames_dispersing, dim = 0, n_eps = 40, n_cores = 1))
t_par <- system.time(cr_par <- crocker(frames_dispersing, dim = 0, n_eps = 40,n_cores = max(1, min(4, n_cores_available))))
cat("Sequential:", t_seq["elapsed"], "s | Parallel:", t_par["elapsed"], "s\n")
