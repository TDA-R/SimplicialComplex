library(Matrix)
library(MASS)
library(gtools)

source("R/Faces.R")
source("R/Boundary.R")
source("R/Betti.R")
source("R/Laplacian.R")

# Example from "Persistent Topological Laplacians—A Survey"
#
#   X = {0,1,2,3, 01,12,23,03}
#   Y = X u {02, 012, 023}
#
# This is the exact example worked through by hand in the paper discussion, so every number below can be checked against that derivation:
#
#   down  (order 01,12,23,03) =
#     [ 2 -1  0  1]
#     [-1  2 -1  0]
#     [ 0 -1  2  1]
#     [ 1  0  1  2]
#
#   upper =
#     [ .5  .5  .5 -.5]
#     [ .5  .5  .5 -.5]
#     [ .5  .5  .5 -.5]
#     [-.5 -.5 -.5  .5]
#
#   sum (= Delta_1^{X,Y}) =
#     [2.5 -.5  .5  .5]
#     [-.5 2.5 -.5 -.5]
#     [ .5 -.5 2.5  .5]
#     [ .5 -.5  .5 2.5]

X_simplices <- list(c(0, 1), c(1, 2), c(2, 3), c(0, 3))
boundary(X_simplices, 1)
Matrix::crossprod(boundary(X_simplices, 1))
Y_simplices <- list(c(0, 1, 2), c(0, 2, 3))
as.matrix(boundary(Y_simplices, 2))

res <- persistent_laplacian(X_simplices, Y_simplices, q = 1)

cat("q-simplex basis (row/column order):", paste(res$labels, collapse = ", "), "\n")

cat("\nDown part  (d_1^X)* d_1^X :\n")
print(round(res$down, 4))

cat("\nUpper part  d_2^{X,Y} (d_2^{X,Y})* :\n")
print(round(res$upper, 4))

cat("\nFull persistent Laplacian  Delta_1^{X,Y} = upper + down :\n")
print(round(res$laplacian, 4))

cat(sprintf("\ndim C_2^{X,Y} = %d  (paper: 1, basis = 023+012)\n", res$dim_upper_chain))

eig <- eigen(res$laplacian, symmetric = TRUE)$values
cat("\nEigenvalues of Delta_1^{X,Y}:\n")
print(round(sort(eig), 4))
cat(sprintf("Zero eigenvalues = persistent Betti number beta_1^{X,Y}: %d\n", sum(abs(eig) < 1e-6)))
cat("(Expected 0: both triangles together fill the whole square, so X's\n", " 1-dimensional loop does not survive into Y.)\n", sep = "")


# X = Y must collapse to the ordinary combinatorial Laplacian,
# whose kernel dimension is just the ordinary Betti number (already computed
# by the package's own betti_number()).

res_XX <- persistent_laplacian(X_simplices, X_simplices, q = 1)
print(res_XX$laplacian)

eig_XX <- eigen(res_XX$laplacian, symmetric = TRUE)$values
cat(sprintf("Zero eigenvalues: %d\n", sum(abs(eig_XX) < 1e-6)))
cat(sprintf("betti_number(X, 1): %d\n", betti_number(X_simplices, 1, tol = 1e-6)))

# Hodge Laplacian
res_hodge <- hodge_laplacian(Y_simplices, k = 1)

cat("k-simplex basis (row/column order):", paste(res_hodge$labels, collapse = ", "), "\n")

cat("\nDown part  B_1^T B_1 :\n")
print(round(res_hodge$down, 4))

cat("\nUp part  B_2 B_2^T :\n")
print(round(res_hodge$up, 4))

cat("\nFull Hodge Laplacian  L_1(Y) = down + up :\n")
print(round(res_hodge$laplacian, 4))

eig_hodge <- eigen(res_hodge$laplacian, symmetric = TRUE)$values
cat("\nEigenvalues of L_1(Y):\n")
print(round(sort(eig_hodge), 4))
cat(sprintf("Zero eigenvalues = beta_1(Y): %d  (betti_number(Y, 1): %d)\n", sum(abs(eig_hodge) < 1e-6), betti_number(Y_simplices, 1, tol = 1e-6)))
cat("(Expected 0: Y is the fully filled-in square, contractible, no loops.)\n")

# Basic graph Laplacian  L = D - A
res_graph <- graph_laplacian(Y_simplices)

cat("\n\n=== Basic graph Laplacian L = D - A  (on Y's 1-skeleton) ===\n")
cat("vertex basis (row/column order):", paste(res_graph$labels, collapse = ", "), "\n")

cat("\nDegree matrix D:\n")
print(res_graph$D)

cat("\nAdjacency matrix A:\n")
print(res_graph$A)

cat("\nGraph Laplacian L = D - A:\n")
print(res_graph$laplacian)

# Sanity check: this must equal hodge_laplacian(Y, 0)$laplacian, since at
# k = 0 the down term vanishes and L_0 = B_1 B_1^T = D - A.
res_hodge0 <- hodge_laplacian(Y_simplices, k = 0)
cat("\nMatches hodge_laplacian(Y_simplices, k = 0)$laplacian? ", isTRUE(all.equal(unname(res_graph$laplacian), unname(res_hodge0$laplacian))), "\n")

eig_graph <- eigen(res_graph$laplacian, symmetric = TRUE)$values
cat("\nEigenvalues of L:\n")
print(round(sort(eig_graph), 4))
cat(sprintf("Zero eigenvalues = beta_0(Y) = number of connected components: %d  (betti_number(Y, 0): %d)\n", sum(abs(eig_graph) < 1e-6), betti_number(Y_simplices, 0, tol = 1e-6)))
