library(testthat)
library(SimplicialComplex)

test_that("faces works correctly", {
  simplices <- list(c(2, 1, 3), c(4, 2), c(5), c(2, 3, 5, 4))
  expect_length(faces(simplices, target_dim=0), 5)
})
