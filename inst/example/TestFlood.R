source("./R/FloodComplex.R")

set.seed(42)
n_pts <- 200000
n_lms <- 250

generate_noisy_torus_points <- function(n, noise_sd = 0.05) {
  R <- 2; r <- 1
  u <- stats::runif(n, 0, 2 * pi)
  # rejection sampling for uniformity on the torus surface
  v <- numeric(n); k <- 0
  while (k < n) {
    cand <- stats::runif(n, 0, 2 * pi)
    acc <- stats::runif(n) < (R + r * cos(cand)) / (R + r)
    take <- min(n - k, sum(acc))
    if (take > 0) { v[(k + 1):(k + take)] <- cand[acc][seq_len(take)]; k <- k + take }
  }
  x <- (R + r * cos(v)) * cos(u)
  y <- (R + r * cos(v)) * sin(u)
  z <- r * sin(v)
  cbind(x, y, z) + matrix(stats::rnorm(3 * n, sd = noise_sd), ncol = 3)
}

pts <- generate_noisy_torus_points(n_pts)

t0 <- Sys.time()
filtration <- build_flood_filtration(pts, landmarks = n_lms, points_per_edge = 15, backend = "cpu")
t1 <- Sys.time()
cat(sprintf("Flood filtration: %d simplices in %.1fs\n", length(filtration), as.numeric(t1 - t0, units = "secs")))

pairs <- flood_persistence(filtration)
t2 <- Sys.time()
cat(sprintf("Persistence: %.1fs\n", as.numeric(t2 - t1, units = "secs")))

# most persistent intervals per dimension
pairs$pers <- pairs$death - pairs$birth
for (d in 0:2) {
  pd <- pairs[pairs$dim == d, ]
  pd <- pd[order(-pd$pers), ]
  cat(sprintf("\nH%d (top 5 of %d):\n", d, nrow(pd)))
  print(head(pd[, c("birth", "death", "pers")], 5), row.names = FALSE)
}

# persistence diagram
fin <- pairs[is.finite(pairs$death), ]
lim <- c(0, max(fin$death) * 1.05)
plot(fin$birth, fin$death, col = fin$dim + 1, pch = 19, cex = 0.7,
     xlim = lim, ylim = lim, xlab = "Birth", ylab = "Death",
     main = sprintf("Flood PH, noisy torus (%d pts, %d landmarks)", n_pts, n_lms))
abline(0, 1, lty = 2, col = "grey")
legend("bottomright", legend = paste0("H", 0:2), col = 1:3, pch = 19)
