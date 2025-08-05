#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

library(gtools)
library(Matrix)

abstract_simplicial_complex <- function(
    simplices, dimension, eps=NULL
) {

  # check
  max_dim <- max(sapply(simplices, length)) -1
  # we need +-1 dimension to calculate boundary and betti number
  # so if the max simplices is 3 point, the max dimension is 2
  if (dimension <= 1 || dimension >= max_dim) {
    stop("dimension must be between 0 and max_dim")
  }
  print(max_dim)

  # simplices
  simplices <- lapply(simplices, function(simplex) sort(simplex))

  # betti numbers
  # betti_num <- betti_number(simplices, dimension, eps)

  euler_characteristic(simplices, dimension, eps)
}

simplices <- list(c(2, 1, 3), c(4, 2), c(5), c(2, 3, 5, 4))
abstract_simplicial_complex(simplices, 2, eps=NULL)

faces <- function(
    simplices, target_dim
    ) {
  faceset <- list()
  facestr_set <- character()

  for (simplex in simplices) {
    node_counts <- length(simplex)

    for (r in node_counts:1) {
      combs <- combinations(n = node_counts, r = r, v = simplex)

      for (i in 1:nrow(combs)) {
        face <- sort(combs[i, ])
        face_str <- paste(face, collapse = "-")

        if (!(face_str %in% facestr_set)) {
          faceset[[length(faceset) + 1]] <- face
          facestr_set <- c(facestr_set, face_str)
        }
      }
    }
  }

  faces_in_dim <- Filter(function(face) length(face) == target_dim + 1, faceset)

  print("Face set")
  # print(facestr_set)
  # print(faceset)
  # print(faces_in_dim)

  return(faces_in_dim)
}

faces(simplices, target_dim=0)
faces(simplices, target_dim=1)
faces(simplices, target_dim=2)
faces(simplices, target_dim=3)

boundary <- function(
    simplices, bound_dim
    ) {
  simplices_dim <- faces(simplices, bound_dim)
  simplices_lower_dim <- faces(simplices, target_dim=(bound_dim-1))

  if (length(simplices_lower_dim) == 0) {

    n <- length(simplices_dim)
    i <- rep(1, n) # row
    j <- 1:n # column
    x <- rep(1, n) # value
    Sparse <- sparseMatrix(i=i, j=j, x=x)

  }
  else {

    # add index to simplices
    keys_dim <- sapply(simplices_dim, function(x) paste(sort(x), collapse= "-"))
    dct_dim <- setNames(seq_along(simplices_dim), keys_dim)

    keys_lower_dim <- sapply(simplices_lower_dim, function(x) paste(sort(x), collapse= "-"))
    dct_lower_dim <- setNames(seq_along(simplices_lower_dim), keys_lower_dim)

    print('dictionary')
    print(dct_dim)
    print('dictionary lower')
    print(dct_lower_dim)

    Sparse <- Matrix(
      0,
      nrow = length(simplices_lower_dim),
      ncol = length(simplices_dim),
      sparse = TRUE
    )

    for (simplex_dim in simplices_dim) {
      for (face in 1:length(simplex_dim)) {
        # remove face in simplex and return (k-1)-simplex
        # \partial_k \sigma = \sum_i (-1)^i [v_0 v_1 \ldots \hat{v}_i \ldots v_k]
        simplex_lower_dim <- c(simplex_dim[-face])
        i <- dct_lower_dim[paste(simplex_lower_dim, collapse = "-")]
        j <- dct_dim[paste(simplex_dim, collapse = "-")]
        print('Lower dimension Simplices')
        print(i)
        print('Simplices')
        print(j)
        # remove odd and even
        Sparse[i, j] <- ifelse(face %% 2 == 1, -1, 1)
      }
    }
    print('Final Sparse')
    print(Sparse)
  }
  return(Sparse)
}

boundary(simplices, 0) # there is no boundary for 0-simplex
boundary(simplices, 1)
boundary(simplices, 2)
boundary(simplices, 3)

betti_number <- function(
    simplices, bound_dim, eps
    ) {
  partial_i <- boundary(simplices, bound_dim)
  partial_i1 <- boundary(simplices, bound_dim + 1)

  if (bound_dim == 0) {
    partial_i_rank <- 0
  }
  else {
    if (is.null(eps)) {
      partial_i_rank <- rankMatrix(partial_i)
    }
    else {
      partial_i_rank <- rankMatrix(partial_i, tol = eps)
    }
  }

  if (is.null(eps)) {
    partial_i1_rank <- rankMatrix(partial_i1)
  }
  else {
    partial_i1_rank <- rankMatrix(partial_i1, tol = eps)
  }

  # \beta_k = \text{dim} (\partial_k) - \text{rank}(\text{Im} \partial_k)) - \text{rank} (\text{Im} \partial_{k+1})
  betti_num <- ncol(partial_i) - partial_i_rank - partial_i1_rank
  print(betti_num)

  return(betti_num)
}

betti_number(simplices, 0, 0.1) # there is no boundary for 0-simplex
betti_number(simplices, 1, 0.1)
betti_number(simplices, 2, 0.1)
betti_number(simplices, 3, 0.1) # this will return error because no 3-simplex to compute

euler_characteristic <- function(
    simplices, eps
    ) {
  max_dim <- max(sapply(simplices, length))
  print(max_dim)
  #　χ = Σ (-1)^k * β_k
  euler_char <- sum(sapply(0:(max_dim), function(bound_dim) {
    (-1)^bound_dim * betti_number(simplices, bound_dim, eps)
  }))

  # χ = Σ (-1)^k * β_k
  # χ = 2: Sphere-like surfaces
  # χ = 1: Disk-like spaces
  # χ = 0: Torus-like or circle-like spaces
  # χ < 0: Surfaces with multiple handles
  return(euler_char)
}

euler_characteristic(simplices, 0.1)

