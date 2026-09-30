library(ggplot2)

source("./R/Persistence.R")
source("./R/DiscreteMorse.R")
source("./R/PlotMorse.R")
source("./R/ComplexUtils.R")

set.seed(42)

nr <- 40L # rows
nc <- 60L # columns
vid <- function(i, j) (i - 1L) * nc + j # 1-based row-major vertex id
coords <- as.matrix(expand.grid(x = seq_len(nc), y = seq_len(nr)))
triangles <- vector("list", 2L * (nr - 1L) * (nc - 1L))
k <- 0L

for (i in seq_len(nr - 1L)) {
  for (j in seq_len(nc - 1L)) {
    v00 <- vid(i, j)
    v10 <- vid(i + 1L, j)
    v01 <- vid(i, j + 1L)
    v11 <- vid(i + 1L, j + 1L)

    k <- k + 1L; triangles[[k]] <- c(v00, v10, v11)
    k <- k + 1L; triangles[[k]] <- c(v00, v01, v11)
  }
}


cat(sprintf("Grid: %d vertices, %d triangles\n", nrow(coords), length(triangles)))

hills <- list(
  list(cx = 15, cy = 10, h = 12, sx = 7, sy = 6),
  list(cx = 30, cy = 25, h = 9,  sx = 6, sy = 7),
  list(cx = 48, cy = 12, h = 10, sx = 7, sy = 6),
  list(cx = 45, cy = 30, h = 7,  sx = 5, sy = 5)
)
nu <- 0.4 # noise
rho <- rep(0, nrow(coords))
for (hl in hills) {
  dx <- coords[, 1] - hl$cx
  dy <- coords[, 2] - hl$cy
  rho <- rho + hl$h * exp(-((dx^2) / (2 * hl$sx^2) + (dy^2) / (2 * hl$sy^2)))
}
rho <- rho + stats::runif(length(rho), 0, nu)

# reconstruct
mr <- morse_recon(triangles, rho, delta = 1.2, vertex_coords = coords)
cat(sprintf("critical (root) vertices: %d\n", length(mr$dmvf$roots)))

plot_morse_recon(mr, show_density = TRUE)

saddle_center <- coords[mr$dmvf$critical_edges$u[1], ]
plot_morse_vpath(mr, center = saddle_center, radius = 10)

#' options(rgl.useNULL = TRUE)
plot_morse_landscape(mr)
#' rgl::rglwidget()

# partial
library(Matrix)
library(gtools)
library(igraph)
library(ggplot2)

source("./R/Faces.R")
source("./R/VRComplex.R")
source("./R/Filtration.R")
pts <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1), ncol = 2)
filist <- build_filtration(pts, method = "VR", eps_max = 1.5, max_dimension = 1)
filtration <- restrict_filtration(filist, 1)
dmvf <- partial_pers_dmvf(filtration, delta = 0)
collect_g(dmvf)

# lower star
triangles <- list(c(1, 2, 3), c(2, 3, 4))
f <- c(0, 1, 1, 2)
filtration <- lower_star_filtration(triangles, f)
filtration
