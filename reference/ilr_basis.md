# Isometric log-ratio (ILR) basis matrix

Builds the default orthonormal ILR basis of Egozcue et al. (2003, eq.
17-18) for a composition with `K` parts: column `i` is the vector
\$\$u_i = \sqrt{i/(i+1)} \\ \[\underbrace{1/i, \dots, 1/i}\_{i}, -1, 0,
\dots, 0\]\$\$ whose `clr`-inverse, `e_i = clr_inv(u_i)`, is the `i`-th
element of the Aitchison-orthonormal simplex basis (eq. 18) contrasting
the (equally weighted) average of the first `i` parts against part
`i + 1` — the same default sequential basis used by `MixSIAR`. `bsimms`
passes this matrix into the generated Stan program as data (`V`), and
recovers source proportions from an unconstrained linear predictor `eta`
(in ILR space) via `softmax(V * eta)` (see
[`ilr_inv()`](https://mattiaghilardi.github.io/bsimms/reference/ilr.md)).

## Usage

``` r
ilr_basis(K)
```

## Arguments

- K:

  Integer, number of parts (sources), `K >= 2`.

## Value

A numeric matrix with `K` rows and `K - 1` columns. Columns are
orthonormal and each column sums to zero, so `V %*% z` always lands on
the zero-sum (clr) hyperplane for any `z`.

## References

Egozcue, J.J., Pawlowsky-Glahn, V., Mateu-Figueras, G., & Barcelo-Vidal,
C. (2003). Isometric logratio transformations for compositional data
analysis. *Mathematical Geology*, 35(3), 279-300.
[doi:10.1023/A:1023818214614](https://doi.org/10.1023/A%3A1023818214614)

## Examples

``` r
V <- ilr_basis(4)
round(crossprod(V), 10) # identity: columns are orthonormal
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    0
#> [3,]    0    0    1
round(colSums(V), 10) # zero: valid clr-constrained basis
#> [1] 0 0 0
```
