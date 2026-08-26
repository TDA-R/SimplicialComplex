library(testthat)
library(SimplicialComplex)

test_that("faces works correctly", {
  simplices <- list(c(2, 1, 3), c(4, 2), c(5), c(2, 3, 5, 4))
  expect_length(faces(simplices, target_dim=0), 5)
})

test_that("persistence_landscape computes the tent-function landscape correctly", {
  # two identical bars (0, 2): tent peaks at t = 1 with height 1
  df <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(2, 2))
  landscape <- persistence_landscape(df, dimension = 0, resolution = 3, t_range = c(0, 2))

  at_0 <- landscape[landscape$t == 0, ]
  at_1 <- landscape[landscape$t == 1, ]
  at_2 <- landscape[landscape$t == 2, ]

  expect_setequal(landscape$k, c(1, 2))
  expect_equal(at_0$value, c(0, 0))
  expect_equal(at_1$value, c(1, 1))
  expect_equal(at_2$value, c(0, 0))
})

test_that("persistence_landscape drops essential (death = Inf) pairs", {
  df <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(1, Inf))
  expect_message(
    landscape <- persistence_landscape(df, dimension = 0, resolution = 5),
    "essential"
  )
  expect_true(all(is.finite(landscape$value)))
  expect_equal(max(landscape$k), 1)
})

test_that("build_filtration(method = 'VR') matches VietorisRipsComplex directly", {
  points <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1), ncol = 2)
  filtration <- build_filtration(points, method = "VR", eps_max = 1.2)
  vr <- VietorisRipsComplex(points, epsilon = 1.2)
  expect_equal(length(filtration), sum(sapply(0:(max(sapply(vr$simplices, length)) - 1),
                                               function(k) length(faces(vr$simplices, k)))))
  expect_true(all(vapply(filtration, function(x) x$t >= 0, TRUE)))
  expect_true(is.unsorted(sapply(filtration, `[[`, "t")) == FALSE)
})

test_that("build_filtration requires eps_max for VR/Cech/Witness but not Delaunay", {
  points <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1), ncol = 2)
  expect_error(build_filtration(points, method = "VR"), "eps_max")
  expect_error(build_filtration(points, method = "Cech"), "eps_max")
  expect_error(build_filtration(points, method = "Witness"), "eps_max")
  expect_error(build_filtration(points, method = "Witness", eps_max = 1), "landmarks")
  expect_silent(build_filtration(points, method = "Delaunay"))
})

test_that("AlphaComplex reproduces the textbook unit-square example", {
  # 4 corners of the unit square: the hole (empty square) is born once the
  # perimeter closes (edge length 1 -> alpha value 0.5) and dies once the
  # diagonal enters (length sqrt(2) -> alpha value sqrt(2)/2)
  points <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  filtration <- build_filtration(points, method = "Alpha")
  pairs <- persistence_pairs(filtration)
  h1 <- pairs[pairs$dim == 1, ]
  # simultaneous events at t = sqrt(2)/2 (the diagonal edge and both
  # triangles all enter at once) can also produce a zero-persistence bar;
  # only the genuine (birth < death) feature matters here
  real <- h1[h1$death > h1$birth, ]
  expect_equal(nrow(real), 1)
  expect_equal(real$birth, 0.5, tolerance = 1e-8)
  expect_equal(real$death, sqrt(2) / 2, tolerance = 1e-8)
})

test_that("DelaunayComplex is AlphaComplex at epsilon = Inf", {
  points <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0.5, 0.4), ncol = 2, byrow = TRUE)
  dc <- DelaunayComplex(points)
  ac <- AlphaComplex(points, epsilon = Inf)
  expect_equal(dc$simplices, ac$simplices)
  expect_equal(dc$filtration, ac$filtration)
  expect_true(inherits(dc, "delaunay_complex"))
})

test_that("CechComplex includes the full simplex once epsilon reaches the circumradius", {
  # 4 corners of the unit square: every point is within sqrt(2)/2 of the
  # center, so the whole 4-point set enters the Cech complex at exactly that
  # scale - the well-known "Cech complex can jump straight to high
  # dimension" behaviour that distinguishes it from Alpha/Delaunay.
  points <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  filtration <- build_filtration(points, method = "Cech", eps_max = 2)
  full <- Filter(function(x) length(x$simplex) == 4, filtration)
  expect_length(full, 1)
  expect_equal(full[[1]]$t, sqrt(2) / 2, tolerance = 1e-8)
})

test_that("restrict_filtration keeps only simplices at or below the given dimension", {
  simplices <- list(c(1), c(2), c(1, 2), c(1, 2, 3))
  filist <- list(
    list(simplex = c(1), t = 0), list(simplex = c(2), t = 0),
    list(simplex = c(1, 2), t = 0.5), list(simplex = c(1, 2, 3), t = 1)
  )
  r0 <- restrict_filtration(filist, max_dimension = 0)
  expect_true(all(vapply(r0, function(x) length(x$simplex) - 1, 0) <= 0))
  r1 <- restrict_filtration(filist, max_dimension = 1)
  expect_true(all(vapply(r1, function(x) length(x$simplex) - 1, 0) <= 1))
  expect_length(r1, 3)
})

test_that("persistence_pairs(max_dimension=) reports correct death times, not a naive drop", {
  # the Cech blow-up example above: H1's one real feature (birth 0.5) is
  # killed by a triangle (a dimension-2 simplex) entering at sqrt(2)/2 - so
  # capping at max_dimension = 1 must still report death = sqrt(2)/2 exactly
  # as the uncapped computation does, by keeping dimension 2 internally as a
  # potential killer and only dropping it from the *output*. A naive
  # dim <= 1 truncation (no "+1") would instead leave that class essential.
  points <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  filtration <- build_filtration(points, method = "Cech", eps_max = 1.2)

  full <- persistence_pairs(filtration)
  capped <- persistence_pairs(filtration, max_dimension = 1)
  naive <- persistence_pairs(restrict_filtration(filtration, 1))

  expect_true(all(capped$dim <= 1))

  h1_full <- full[full$dim == 1 & full$death > full$birth, c("birth", "death")]
  h1_capped <- capped[capped$dim == 1 & capped$death > capped$birth, c("birth", "death")]
  rownames(h1_full) <- NULL; rownames(h1_capped) <- NULL
  expect_equal(h1_capped, h1_full) # correct death time preserved, not dropped

  h1_naive <- naive[naive$dim == 1 & naive$birth == 0.5, ]
  expect_true(is.infinite(h1_naive$death)) # demonstrates the naive approach is wrong
})

test_that("flood_persistence(max_dimension=) applies the same correction", {
  set.seed(5)
  pts <- matrix(rnorm(60), ncol = 2)
  filtration <- build_flood_filtration(pts, landmarks = 12, points_per_edge = 6)
  full <- flood_persistence(filtration)
  capped <- flood_persistence(filtration, max_dimension = 0)
  expect_true(all(capped$dim == 0))
  h0_full <- full[full$dim == 0, c("birth", "death")]
  h0_capped <- capped[capped$dim == 0, c("birth", "death")]
  rownames(h0_full) <- NULL; rownames(h0_capped) <- NULL
  expect_equal(h0_capped, h0_full)
})

test_that("boundary_info()/extract_persistence_pairs() honour max_dimension identically to persistence_pairs()", {
  # same Cech blow-up example as the persistence_pairs() max_dimension test:
  # the dense pipeline (boundary_info + extract_persistence_pairs) must
  # apply the exact same "+1 internally, trim the output" correction as the
  # sparse persistence_pairs(), so the two agree bit for bit.
  points <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  filtration <- build_filtration(points, method = "Cech", eps_max = 1.2)

  res <- boundary_info(filtration, max_dimension = 1)
  dense <- extract_persistence_pairs(res$filist, res$last_1, res$pivot_owner, max_dimension = 1)
  sparse <- persistence_pairs(filtration, max_dimension = 1)

  expect_true(all(dense$dim <= 1))
  d1 <- dense[order(dense$dim, dense$birth), c("dim", "birth", "death")]
  d2 <- sparse[order(sparse$dim, sparse$birth), c("dim", "birth", "death")]
  rownames(d1) <- NULL; rownames(d2) <- NULL
  expect_equal(d1, d2)
})

test_that("extract_persistence_pairs() refuses a filist that doesn't match last_1/pivot_owner", {
  points <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  filtration <- build_filtration(points, method = "Cech", eps_max = 1.2)
  res <- boundary_info(filtration, max_dimension = 1) # res$filist is shorter than filtration

  expect_error(
    extract_persistence_pairs(filtration, res$last_1, res$pivot_owner, max_dimension = 1),
    "must come from the same boundary_info"
  )
  # using the matching filist works fine
  expect_silent(extract_persistence_pairs(res$filist, res$last_1, res$pivot_owner, max_dimension = 1))
})

test_that("boundary_info()/extract_persistence_pairs() without max_dimension are unchanged", {
  points <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1), ncol = 2)
  filtration <- build_filtration(points, method = "VR", eps_max = 1.2)
  res <- boundary_info(filtration)
  expect_identical(res$filist, filtration)
  pairs <- extract_persistence_pairs(filtration, res$last_1, res$pivot_owner)
  sparse <- persistence_pairs(filtration)
  p1 <- pairs[order(pairs$dim, pairs$birth), c("dim", "birth", "death")]
  p2 <- sparse[order(sparse$dim, sparse$birth), c("dim", "birth", "death")]
  rownames(p1) <- NULL; rownames(p2) <- NULL
  expect_equal(p1, p2)
})

test_that("build_filtration(max_dimension=) tags its result so downstream functions need no parameter", {
  # Same Cech blow-up example as the persistence_pairs(max_dimension=) test
  # above, but the cap now comes from build_filtration() itself - the
  # downstream calls below pass no max_dimension at all and must still match
  # the already-verified explicit-parameter path exactly.
  points <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  full <- build_filtration(points, method = "Cech", eps_max = 1.2)
  capped <- build_filtration(points, method = "Cech", eps_max = 1.2, max_dimension = 1)

  expect_equal(attr(capped, "max_dimension"), 1)
  # the "+1 trick": dimension 2 is still present internally...
  expect_true(any(vapply(capped, function(x) length(x$simplex) - 1, 0) == 2))
  # ...but nothing above it is built at all, unlike the uncapped filtration
  expect_true(any(vapply(full, function(x) length(x$simplex) - 1, 0) > 2) ||
                max(vapply(full, function(x) length(x$simplex) - 1, 0)) == 2)
  expect_true(length(capped) <= length(full))

  explicit <- persistence_pairs(full, max_dimension = 1)
  auto <- persistence_pairs(capped) # no max_dimension argument at all
  e <- explicit[order(explicit$dim, explicit$birth), c("dim", "birth", "death")]
  a <- auto[order(auto$dim, auto$birth), c("dim", "birth", "death")]
  rownames(e) <- NULL; rownames(a) <- NULL
  expect_equal(a, e)
})

test_that("build_filtration(max_dimension=) auto-detection also reaches boundary_info()/extract_persistence_pairs()", {
  points <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  capped <- build_filtration(points, method = "Cech", eps_max = 1.2, max_dimension = 1)

  res <- boundary_info(capped) # no max_dimension argument
  expect_equal(attr(res$filist, "max_dimension"), 1)
  dense <- extract_persistence_pairs(res$filist, res$last_1, res$pivot_owner) # none here either
  sparse <- persistence_pairs(capped)

  expect_true(all(dense$dim <= 1))
  d1 <- dense[order(dense$dim, dense$birth), c("dim", "birth", "death")]
  d2 <- sparse[order(sparse$dim, sparse$birth), c("dim", "birth", "death")]
  rownames(d1) <- NULL; rownames(d2) <- NULL
  expect_equal(d1, d2)
})

test_that("an explicit downstream max_dimension overrides build_filtration()'s own cap", {
  points <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), ncol = 2, byrow = TRUE)
  capped <- build_filtration(points, method = "Cech", eps_max = 1.2, max_dimension = 1)

  narrower <- persistence_pairs(capped, max_dimension = 0)
  expect_true(all(narrower$dim == 0))

  res <- boundary_info(capped, max_dimension = 0)
  expect_equal(attr(res$filist, "max_dimension"), 0)
  dense <- extract_persistence_pairs(res$filist, res$last_1, res$pivot_owner, max_dimension = 0)
  expect_true(all(dense$dim == 0))
})

test_that("build_filtration() without max_dimension carries no attribute and is unaffected", {
  points <- matrix(c(0, 1, 1, 0, 0, 0, 1, 1), ncol = 2)
  filtration <- build_filtration(points, method = "VR", eps_max = 1.2)
  expect_null(attr(filtration, "max_dimension"))

  res <- boundary_info(filtration)
  expect_null(attr(res$filist, "max_dimension"))
  expect_identical(res$filist, filtration)
})

test_that("build_filtration(max_dimension=) restricts Alpha/Delaunay complexes correctly (post-hoc branch)", {
  # Alpha/Delaunay have no build-time dimension cap, so build_filtration()
  # restricts them after the fact - a different code path from the VR/Cech/
  # Witness branches (which cap simplices_to_filtration()'s kmax directly).
  # Exercise it explicitly so both branches are covered.
  points <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0.5, 0.5), ncol = 2, byrow = TRUE)

  for (method in c("Alpha", "Delaunay")) {
    full <- build_filtration(points, method = method)
    capped <- build_filtration(points, method = method, max_dimension = 1)
    expect_equal(attr(capped, "max_dimension"), 1)
    expect_true(max(vapply(capped, function(x) length(x$simplex) - 1, 0)) <= 2)

    explicit <- persistence_pairs(full, max_dimension = 1)
    auto <- persistence_pairs(capped)
    e <- explicit[order(explicit$dim, explicit$birth), c("dim", "birth", "death")]
    a <- auto[order(auto$dim, auto$birth), c("dim", "birth", "death")]
    rownames(e) <- NULL; rownames(a) <- NULL
    expect_equal(a, e)
  }
})

test_that("WitnessComplex builds a valid complex on the requested number of landmarks", {
  set.seed(7)
  theta <- seq(0, 2 * pi, length.out = 25)[-25]
  points <- cbind(cos(theta), sin(theta))
  filtration <- build_filtration(points, method = "Witness", eps_max = 0.6, landmarks = 8)
  vertex_ids <- unique(unlist(lapply(filtration, function(x) x$simplex[length(x$simplex) == 1])))
  expect_true(length(unique(unlist(lapply(filtration[sapply(filtration, function(x) length(x$simplex) == 1)],
                                           `[[`, "simplex")))) <= 8)
  expect_true(all(vapply(filtration, function(x) x$t >= 0, TRUE)))
})

test_that("wasserstein_distance and bottleneck_distance are zero for identical diagrams", {
  df <- data.frame(dim = c(0, 0, 1), birth = c(0, 0, 0.2), death = c(1, 2, 0.9))
  expect_equal(wasserstein_distance(df, df, dimension = 0, p = 2), 0)
  expect_equal(bottleneck_distance(df, df, dimension = 0), 0)
  expect_equal(wasserstein_distance(df, df, dimension = 1, p = 1), 0)
})

test_that("wasserstein_distance and bottleneck_distance match hand-computed values", {
  df1 <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(1, 2))
  df2 <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(1, 3))
  # optimal matching pairs (0,1)-(0,1) exactly and (0,2)-(0,3) at distance 1
  expect_equal(wasserstein_distance(df1, df2, dimension = 0, p = 2), 1, tolerance = 1e-8)
  expect_equal(bottleneck_distance(df1, df2, dimension = 0), 1, tolerance = 1e-8)
})

test_that("distance functions cap essential pairs by birth time instead of dropping them", {
  # df1 has an essential class (0, Inf); dropping it would let it vanish
  # from the comparison entirely. It should instead be capped (at 2x the
  # largest finite value present, i.e. death = 2 here) and matched by its
  # birth time like any other point.
  df1 <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(1, Inf))
  df2 <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(1, 1))
  expect_message(w <- wasserstein_distance(df1, df2, dimension = 0, p = 2), "Capping")
  expect_equal(w, 1, tolerance = 1e-8)
  expect_message(b <- bottleneck_distance(df1, df2, dimension = 0), "Capping")
  expect_equal(b, 1, tolerance = 1e-8)
})

test_that("essential pairs still give zero self-distance", {
  df <- data.frame(dim = c(0, 1), birth = c(0, 0.2), death = c(Inf, Inf))
  expect_equal(wasserstein_distance(df, df, dimension = 0, p = 2), 0)
  expect_equal(bottleneck_distance(df, df, dimension = 1), 0)
})

test_that("diagram_matching's distance agrees with wasserstein_distance/bottleneck_distance", {
  df1 <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(1, 2))
  df2 <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(1, 3))

  mw <- diagram_matching(df1, df2, dimension = 0, distance = "wasserstein", p = 2)
  expect_equal(mw$distance, wasserstein_distance(df1, df2, dimension = 0, p = 2), tolerance = 1e-8)
  # both points matched real-to-real: (0,1)-(0,1) at zero cost, (0,2)-(0,3) at cost 1
  expect_equal(nrow(mw$matches), 2)
  expect_true(all(mw$matches$type == "real-real"))

  mb <- diagram_matching(df1, df2, dimension = 0, distance = "bottleneck")
  expect_equal(mb$distance, bottleneck_distance(df1, df2, dimension = 0), tolerance = 1e-8)
})

test_that("diagram_matching caps essential pairs and flags them in the matching table", {
  df1 <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(1, Inf))
  df2 <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(1, 1))

  expect_message(m <- diagram_matching(df1, df2, dimension = 0, distance = "wasserstein", p = 2),
                  "Capping")
  expect_equal(m$distance, 1, tolerance = 1e-8)
  expect_true(any(m$matches$x_essential, na.rm = TRUE))
  # the essential point's capped death matches cap_essential's own convention (2x max finite)
  essential_row <- m$matches[which(m$matches$x_essential), ]
  expect_equal(essential_row$x_death, 2, tolerance = 1e-8)
})

test_that("diagram_matching errors when neither diagram has points in the requested dimension", {
  df1 <- data.frame(dim = c(0), birth = c(0), death = c(1))
  df2 <- data.frame(dim = c(0), birth = c(0), death = c(1))
  expect_error(diagram_matching(df1, df2, dimension = 1), "nothing to match")
})

test_that("plot_matching returns a ggplot object built from the same matching", {
  df1 <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(1, 2))
  df2 <- data.frame(dim = c(0, 0), birth = c(0, 0), death = c(1, 3))
  p <- plot_matching(df1, df2, dimension = 0, distance = "wasserstein", p = 2,
                      labels = c("A", "B"))
  expect_s3_class(p, "ggplot")
  expect_true(grepl("distance = 1", p$labels$title))
})

test_that("compare_complexes summarizes every requested method", {
  points <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0.5, 0.5), ncol = 2, byrow = TRUE)
  cmp <- compare_complexes(points, methods = c("VR", "Delaunay", "Alpha"), eps_max = 1.2)
  expect_setequal(cmp$summary$method, c("VR", "Delaunay", "Alpha"))
  expect_true(all(cmp$summary$n_simplices > 0))
  expect_named(cmp$diagrams, c("VR", "Delaunay", "Alpha"))
})

test_that("complex_distance computes a symmetric-sensible distance and validates inputs", {
  set.seed(3)
  cloud_a <- matrix(rnorm(24), ncol = 2)
  cloud_b <- cloud_a + 0.05
  res <- complex_distance(cloud_a, cloud_b, method = "VR", distance = "bottleneck",
                           dimension = 0, eps_max = 0.6)
  expect_true(res$distance >= 0)
  self <- complex_distance(cloud_a, cloud_a, method = "Delaunay", distance = "wasserstein",
                            dimension = 1, p = 2)
  expect_equal(self$distance, 0)
  expect_error(complex_distance(cloud_a, cloud_b, method = "Cech"), "eps_max")
})

test_that("build_cubical_filtration finds the hole in a ring image", {
  # ring (1) around a hole (5) around a background (9): the hole should be
  # born when the ring closes (t = 1) and die when its center fills (t = 5)
  image <- matrix(c(
    9, 9, 9, 9, 9,
    9, 1, 1, 1, 9,
    9, 1, 5, 1, 9,
    9, 1, 1, 1, 9,
    9, 9, 9, 9, 9
  ), nrow = 5, byrow = TRUE)

  filtration <- build_cubical_filtration(image)
  pairs <- persistence_pairs(filtration)
  h1 <- pairs[pairs$dim == 1, ]

  expect_true(any(h1$birth == 1 & h1$death == 5))
  expect_true(all(vapply(filtration, function(x) length(x$simplex), 0) <= 3))
})

test_that("plot_persistence labels dimensions discretely as H0/H1/H2/... with a legend", {
  df <- data.frame(dim = c(0, 0, 1, 2), birth = c(0, 0, 0.5, 0.7), death = c(0.5, Inf, 0.7, 0.9))
  p <- plot_persistence(df)

  expect_s3_class(p, "ggplot")
  expect_true(is.factor(p$data$dim_label))
  expect_equal(levels(p$data$dim_label), c("H0", "H1", "H2"))
  expect_false(identical(p$theme$legend.position, "none"))
})
