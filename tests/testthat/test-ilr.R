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

test_that("ilr_basis errors on an invalid K", {
  expect_snapshot(error = TRUE, ilr_basis(K = 1))
  expect_snapshot(error = TRUE, ilr_basis(K = NA))
})

test_that("clr / clr_inv round-trip on a single composition", {
  p <- c(0.5, 0.3, 0.2)
  z <- clr(p)
  expect_length(z, 3)
  expect_equal(sum(z), 0, tolerance = 1e-10)
  expect_equal(clr_inv(z), p, tolerance = 1e-8)
})

test_that("clr / clr_inv round-trip on a matrix of compositions", {
  set.seed(1)
  P <- matrix(runif(30), nrow = 10)
  P <- P / rowSums(P)
  Z <- clr(P)
  expect_equal(dim(Z), c(10, 3))
  expect_equal(unname(rowSums(Z)), rep(0, 10), tolerance = 1e-10)
  expect_equal(clr_inv(Z), P, tolerance = 1e-8)
})

test_that("clr errors on non-positive parts, vector and matrix", {
  expect_snapshot(error = TRUE, clr(c(0.5, -0.1, 0.6)))
  expect_snapshot(error = TRUE, clr(matrix(c(0.5, 0, 0.3, 0.2), nrow = 2)))
})

test_that("ilr errors on non-positive parts, vector and matrix (default basis)", {
  expect_snapshot(error = TRUE, ilr(c(0.5, -0.1, 0.6)))
  expect_snapshot(error = TRUE, ilr(matrix(c(0.5, 0, 0.3, 0.2), nrow = 2)))
})

test_that("ilr with an explicit V matches the default sequential basis", {
  p <- c(0.5, 0.3, 0.2)
  V <- ilr_basis(K = 3)
  expect_equal(ilr(x = p, V = V), ilr(x = p), tolerance = 1e-10)
  expect_equal(ilr_inv(z = ilr(x = p, V = V), V = V), p, tolerance = 1e-8)
})

test_that("ilr with an explicit V round-trips on a matrix of compositions", {
  set.seed(1)
  P <- matrix(runif(30), nrow = 10)
  P <- P / rowSums(P)
  V <- ilr_basis(K = 3)
  Z <- ilr(x = P, V = V)
  expect_equal(dim(Z), c(10, 2))
  expect_equal(ilr_inv(z = Z, V = V), P, tolerance = 1e-8)
})
