# Centred log-ratio transform and its inverse

`clr()` maps compositions to the (mean-centred) log scale; `clr_inv()`
(a softmax) maps back to the simplex.

## Usage

``` r
clr(x)

clr_inv(y)
```

## Arguments

- x:

  A vector or matrix of strictly positive compositions. If a matrix,
  rows are compositions.

- y:

  A vector or matrix on the clr scale (rows are clr vectors if a
  matrix).

## Value

`clr(x)`, same shape as `x`. `clr_inv(y)` returns a composition (or
matrix of compositions) on the simplex.

## Examples

``` r
p <- c(0.5, 0.3, 0.2)
z <- clr(p)
clr_inv(z) # back to p
#> [1] 0.5 0.3 0.2
```
