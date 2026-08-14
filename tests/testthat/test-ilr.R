test_that("ilr_basis is orthonormal and zero-sum", {
  for (K in 2:6) {
    V <- ilr_basis(K = K)
    expect_equal(dim(V), c(K, K - 1))
    expect_equal(unname(crossprod(V)), diag(K - 1), tolerance = 1e-10)
    expect_equal(unname(colSums(V)), rep(0, K - 1), tolerance = 1e-10)
  }
})

test_that("ilr / ilr_inv round-trip on a single composition", {
  p <- c(0.5, 0.3, 0.2)
  z <- ilr(x = p)
  expect_length(z, 2)
  expect_equal(ilr_inv(z = z), p, tolerance = 1e-8)
})

test_that("ilr / ilr_inv round-trip on a matrix of compositions", {
  set.seed(1)
  P <- matrix(runif(30), nrow = 10)
  P <- P / rowSums(P)
  Z <- ilr(x = P)
  expect_equal(dim(Z), c(10, 2))
  expect_equal(ilr_inv(z = Z), P, tolerance = 1e-8)
})

test_that("ilr is invariant to the overall scale of x", {
  p <- c(3, 2, 1)
  expect_equal(ilr(x = p), ilr(x = p / sum(p)))
})
