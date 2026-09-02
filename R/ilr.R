#' Isometric log-ratio (ILR) basis matrix
#'
#' Builds the default orthonormal ILR basis of Egozcue et al. (2003,
#' eq. 17-18) for a composition with `K` parts: column `i` is the vector
#' \deqn{u_i = \sqrt{i/(i+1)} \, [\underbrace{1/i, \dots, 1/i}_{i}, -1, 0,
#' \dots, 0]}
#' whose `clr`-inverse, `e_i = clr_inv(u_i)`, is the `i`-th element of the
#' Aitchison-orthonormal simplex basis (eq. 18) contrasting the (equally
#' weighted) average of the first `i` parts against part `i + 1` — the
#' same default sequential basis used by `MixSIAR`. `bsimms` passes this
#' matrix into the generated Stan program as data (`V`); the Stan model
#' derives `E`, the simplex-domain basis (`E[, d] = softmax(V[, d])`), in
#' `transformed data`, and recovers source proportions from an
#' unconstrained linear predictor `eta` (in ILR space) via the
#' `inverse_ilr()` Stan function, implementing eq. 24 (see [ilr_inv()]).
#'
#' @param K Integer, number of parts (sources), `K >= 2`.
#' @return A numeric matrix with `K` rows and `K - 1` columns. Columns are
#'   orthonormal and each column sums to zero, so `V %*% z` always lands on
#'   the zero-sum (clr) hyperplane for any `z`.
#' @references Egozcue, J.J., Pawlowsky-Glahn, V., Mateu-Figueras, G., &
#'   Barcelo-Vidal, C. (2003). Isometric logratio transformations for
#'   compositional data analysis. *Mathematical Geology*, 35(3), 279-300.
#'   \doi{10.1023/A:1023818214614}
#' @export
#' @examples
#' V <- ilr_basis(4)
#' round(crossprod(V), 10)   # identity: columns are orthonormal
#' round(colSums(V), 10)     # zero: valid clr-constrained basis
ilr_basis <- function(K) {
  K <- as.integer(K)
  if (is.na(K) || K < 2) {
    cli::cli_abort("{.arg K} must be an integer >= 2.", call = NULL)
  }
  V <- matrix(0, nrow = K, ncol = K - 1)
  for (j in seq_len(K - 1)) {
    cst <- sqrt(j / (j + 1))
    V[seq_len(j), j] <- cst / j
    V[j + 1, j] <- -cst
  }
  V
}

#' Centred log-ratio transform and its inverse
#'
#' `clr()` maps compositions to the (mean-centred) log scale; `clr_inv()`
#' (a softmax) maps back to the simplex.
#'
#' @param x A vector or matrix of strictly positive compositions. If a
#'   matrix, rows are compositions.
#' @param y A vector or matrix on the clr scale (rows are clr vectors if a
#'   matrix).
#' @return `clr(x)`, same shape as `x`. `clr_inv(y)` returns a composition
#'   (or matrix of compositions) on the simplex.
#' @export
#' @examples
#' p <- c(0.5, 0.3, 0.2)
#' z <- clr(p)
#' clr_inv(z)  # back to p
clr <- function(x) {
  if (is.matrix(x)) {
    if (any(x <= 0)) cli::cli_abort("All parts of {.arg x} must be strictly positive.", call = NULL)
    lx <- log(x)
    lx - rowMeans(lx)
  } else {
    if (any(x <= 0)) cli::cli_abort("All parts of {.arg x} must be strictly positive.", call = NULL)
    lx <- log(x)
    lx - mean(lx)
  }
}

#' @rdname clr
#' @export
clr_inv <- function(y) {
  if (is.matrix(y)) {
    ey <- exp(y - apply(y, 1, max))
    ey / rowSums(ey)
  } else {
    ey <- exp(y - max(y))
    ey / sum(ey)
  }
}

#' Forward and inverse ILR transforms
#'
#' `ilr()` maps compositions to isometric log-ratio (ILR) coordinates;
#' `ilr_inv()` maps back to the simplex.
#'
#' For the default sequential basis (`V = NULL`), `ilr()` computes each
#' coordinate directly from ratios of geometric means, following Egozcue
#' et al. (2003, eq. 25) — the same formula used by `MixSIAR`:
#' \deqn{y_i = \sqrt{i/(i+1)} \, \ln\!\left(\frac{g(x_1, \dots,
#' x_i)}{x_{i+1}}\right), \quad i = 1, \dots, K - 1,}
#' where \eqn{g(\cdot)} denotes the geometric mean. For a non-default
#' basis, the general definition (eq. 23), `ilr(x) = crossprod(V, clr(x))`,
#' is used instead.
#'
#' `ilr_inv()` implements the inverse transformation (eq. 24), \eqn{x =
#' \bigoplus_{i} (y_i \otimes e_i)}: each ILR coordinate `z[d]`
#' perturbation-scales its simplex-domain basis element `e_d =
#' clr_inv(V[, d])` (a power operation followed by closure/normalisation),
#' and the results are combined by repeated perturbation (elementwise
#' product followed by closure). This is the same construction used by
#' `MixSIAR`'s JAGS implementation, and is what `bsimms`'s generated Stan
#' code uses internally (the `inverse_ilr()` function in the `functions`
#' block of `make_stancode()`'s output) to map the ILR-scale linear
#' predictor back onto the source simplex.
#'
#' @param x A vector of length `K`, or an `N x K` matrix, of strictly
#'   positive parts (a composition, or set of compositions, on the `K`-part
#'   simplex). Need not sum to 1: `ilr()` is invariant to the overall scale
#'   of `x` (only the relative proportions matter), so an unnormalised
#'   vector of positive weights works just as well as a closed composition.
#' @param z A vector of length `K - 1`, or an `N x (K - 1)` matrix of ILR
#'   coordinates.
#' @param V Optional ILR basis from [ilr_basis()]. If `NULL` (default):
#'   for `ilr()`, the default sequential basis (eq. 18) is used, via the
#'   closed-form eq. 25; for `ilr_inv()`, computed automatically from
#'   `ncol(z) + 1` / `length(z) + 1`.
#' @return `ilr()` returns ILR coordinates: a vector of length `K - 1`, or
#'   an `N x (K - 1)` matrix. `ilr_inv()` returns proportions on the
#'   `K`-part simplex: a vector of length `K`, or an `N x K` matrix.
#' @references Egozcue, J.J., Pawlowsky-Glahn, V., Mateu-Figueras, G., &
#'   Barcelo-Vidal, C. (2003). Isometric logratio transformations for
#'   compositional data analysis. *Mathematical Geology*, 35(3), 279-300.
#'   \doi{10.1023/A:1023818214614}
#' @export
#' @examples
#' p <- c(0.5, 0.3, 0.2)
#' z <- ilr(p)
#' ilr_inv(z)  # back to p
ilr <- function(x, V = NULL) {
  if (!is.null(V)) {
    cx <- clr(x)
    return(if (is.matrix(x)) cx %*% V else as.numeric(crossprod(V, cx)))
  }
  if (is.matrix(x)) {
    K <- ncol(x)
    t(vapply(seq_len(nrow(x)), function(i) ilr_default_vec(x[i, ]), numeric(K - 1L)))
  } else {
    ilr_default_vec(x)
  }
}

#' Forward ILR transform (default sequential basis, eq. 25) for a single
#' composition. The per-row worker behind [ilr()]'s `V = NULL` path.
#'
#' @param x Numeric vector of length `K`, strictly positive parts.
#' @return Numeric vector of length `K - 1`: ILR coordinates of `x`.
#' @noRd
ilr_default_vec <- function(x) {
  if (any(x <= 0)) cli::cli_abort("All parts of {.arg x} must be strictly positive.", call = NULL)
  K <- length(x)
  z <- numeric(K - 1L)
  for (i in seq_len(K - 1L)) {
    g_i <- exp(mean(log(x[seq_len(i)])))
    z[i] <- sqrt(i / (i + 1)) * log(g_i / x[i + 1])
  }
  z
}

#' @rdname ilr
#' @importFrom rlang %||%
#' @export
ilr_inv <- function(z, V = NULL) {
  D <- if (is.matrix(z)) ncol(z) else length(z)
  K <- D + 1L
  V <- V %||% ilr_basis(K)
  E <- vapply(seq_len(D), function(d) clr_inv(V[, d]), numeric(K))   # e_1, ..., e_D (eq. 18)
  if (is.matrix(z)) {
    t(vapply(seq_len(nrow(z)), function(i) ilr_inv_vec(z[i, ], E), numeric(K)))
  } else {
    ilr_inv_vec(z, E)
  }
}

#' Inverse ILR transform (eq. 24) for a single ILR coordinate vector, given
#' a precomputed simplex-domain basis. The per-row worker behind
#' [ilr_inv()].
#'
#' @param z Numeric vector of length `D` (`K - 1`), ILR coordinates.
#' @param E Numeric `K x D` simplex-domain basis matrix (columns
#'   `e_1, ..., e_D`, e.g. as built from `V` in [ilr_inv()]).
#' @return Numeric vector of length `K`: a composition on the simplex.
#' @noRd
ilr_inv_vec <- function(z, E) {
  K <- nrow(E)
  D <- ncol(E)
  cross <- matrix(0, K, D)
  for (d in seq_len(D)) {
    powered <- E[, d]^z[d]
    cross[, d] <- powered / sum(powered)
  }
  p <- cross[, 1]
  if (D > 1) {
    for (d in 2:D) p <- p * cross[, d]
  }
  p / sum(p)
}
