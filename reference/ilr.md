# Forward and inverse ILR transforms

`ilr()` maps compositions to isometric log-ratio (ILR) coordinates;
`ilr_inv()` maps back to the simplex.

## Usage

``` r
ilr(x, V = NULL)

ilr_inv(z, V = NULL)
```

## Arguments

- x:

  A vector of length `K`, or an `N x K` matrix, of strictly positive
  parts (a composition, or set of compositions, on the `K`-part
  simplex). Need not sum to 1: `ilr()` is invariant to the overall scale
  of `x` (only the relative proportions matter), so an unnormalised
  vector of positive weights works just as well as a closed composition.

- V:

  Optional ILR basis from
  [`ilr_basis()`](https://mattiaghilardi.github.io/bsimms/reference/ilr_basis.md).
  If `NULL` (default): for `ilr()`, the default sequential basis
  (eq. 18) is used, via the closed-form eq. 25; for `ilr_inv()`,
  computed automatically from `ncol(z) + 1` / `length(z) + 1`.

- z:

  A vector of length `K - 1`, or an `N x (K - 1)` matrix of ILR
  coordinates.

## Value

`ilr()` returns ILR coordinates: a vector of length `K - 1`, or an
`N x (K - 1)` matrix. `ilr_inv()` returns proportions on the `K`-part
simplex: a vector of length `K`, or an `N x K` matrix.

## Details

For the default sequential basis (`V = NULL`), `ilr()` computes each
coordinate directly from ratios of geometric means, following Egozcue et
al. (2003, eq. 25) — the same formula used by `MixSIAR`: \$\$y_i =
\sqrt{i/(i+1)} \\ \ln\\\left(\frac{g(x_1, \dots, x_i)}{x\_{i+1}}\right),
\quad i = 1, \dots, K - 1,\$\$ where \\g(\cdot)\\ denotes the geometric
mean. For a non-default basis, the general definition (eq. 23),
`ilr(x) = crossprod(V, clr(x))`, is used instead.

`ilr_inv()` implements the inverse transformation (eq. 24), \\x =
\bigoplus\_{i} (y_i \otimes e_i)\\, computed in its closed form
`clr_inv(V %*% z)`: a single, numerically stable softmax, algebraically
equivalent to the perturbation-and-closure construction but without its
overflow/underflow risk for large `z`. This is what `bsimms`'s generated
Stan code uses internally (via the built-in `softmax()` function in
[`make_stancode()`](https://mattiaghilardi.github.io/bsimms/reference/make_stancode.md)'s
output) to map the ILR-scale linear predictor back onto the source
simplex.

## References

Egozcue, J.J., Pawlowsky-Glahn, V., Mateu-Figueras, G., & Barcelo-Vidal,
C. (2003). Isometric logratio transformations for compositional data
analysis. *Mathematical Geology*, 35(3), 279-300.
[doi:10.1023/A:1023818214614](https://doi.org/10.1023/A%3A1023818214614)

## Examples

``` r
p <- c(0.5, 0.3, 0.2)
z <- ilr(p)
ilr_inv(z) # back to p
#> [1] 0.5 0.3 0.2
```
