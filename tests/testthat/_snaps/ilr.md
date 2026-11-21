# ilr_basis errors on an invalid K

    Code
      ilr_basis(K = 1)
    Condition
      Error:
      ! `K` must be an integer >= 2.

---

    Code
      ilr_basis(K = NA)
    Condition
      Error:
      ! `K` must be an integer >= 2.

# clr errors on non-positive parts, vector and matrix

    Code
      clr(c(0.5, -0.1, 0.6))
    Condition
      Error:
      ! All parts of `x` must be strictly positive.

---

    Code
      clr(matrix(c(0.5, 0, 0.3, 0.2), nrow = 2))
    Condition
      Error:
      ! All parts of `x` must be strictly positive.

# ilr errors on non-positive parts, vector and matrix (default basis)

    Code
      ilr(c(0.5, -0.1, 0.6))
    Condition
      Error:
      ! All parts of `x` must be strictly positive.

---

    Code
      ilr(matrix(c(0.5, 0, 0.3, 0.2), nrow = 2))
    Condition
      Error:
      ! All parts of `x` must be strictly positive.

