# Compare `bsimms` models by approximate leave-one-out cross-validation

Computes [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html) or
[`loo::waic()`](https://mc-stan.org/loo/reference/waic.html) for each
model (a `loo`/`waic` object passed directly, e.g. one already computed
via `loo(fit)`, is used as-is; for a `bsimms_fit`, a criterion already
cached via
[`add_criterion()`](https://mattiaghilardi.github.io/bsimms/reference/add_criterion.md)
is reused instead of recomputed) and ranks them with
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html).

## Usage

``` r
# S3 method for class 'bsimms_fit'
loo_compare(x, ..., criterion = c("loo", "waic"), model_names = NULL)
```

## Arguments

- x:

  A `bsimms_fit` object. Must be a `bsimms_fit` (not a precomputed
  `loo`/`waic` object) for this method to be dispatched to; a
  precomputed `loo`/`waic` object can only appear in `...`.

- ...:

  Additional `bsimms_fit` objects and/or `loo`/`waic` objects to compare
  against `x`.

- criterion:

  `"loo"` (default) or `"waic"`; which criterion to use for any
  `bsimms_fit` objects passed in (cached via
  [`add_criterion()`](https://mattiaghilardi.github.io/bsimms/reference/add_criterion.md)
  if present, computed on the fly otherwise). Ignored for arguments that
  are already `loo`/`waic` objects.

- model_names:

  Optional character vector naming the compared models, in the order
  `x`, `...`. Defaults to the deparsed argument expressions.

## Value

The result of
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html):
model names are given in its `model` column, ranked by expected log
pointwise predictive density (best model first).

## Examples

``` r
# \donttest{
sim <- simulate_bsimms_data(
  ~Sex,
  n_mixture_obs = 10,
  source_names = c("Beaver", "Deer", "Hare"),
  isotope_names = c("d13C", "d15N"),
  source_means_sds = TRUE,
  n_levels = list(Sex = 2),
  seed = 1
)
fit1 <- bsimm(
  ~1,
  mixture_data = sim$mixture_data,
  source_data = sim$source_data,
  tdf_data = sim$tdf_data,
  isotope_names = sim$isotope_names,
  source_means_sds = sim$source_means_sds,
  tdf_means_sds = sim$tdf_means_sds,
  conc_dep = sim$conc_dep,
  error_structure = sim$error_structure,
  source_col = sim$source_col,
  chains = 2,
  iter_warmup = 500,
  iter_sampling = 500
)
#> Running MCMC with 2 sequential chains...
#> 
#> Chain 1 Iteration:   1 / 1000 [  0%]  (Warmup) 
#> Chain 1 Iteration: 100 / 1000 [ 10%]  (Warmup) 
#> Chain 1 Iteration: 200 / 1000 [ 20%]  (Warmup) 
#> Chain 1 Iteration: 300 / 1000 [ 30%]  (Warmup) 
#> Chain 1 Iteration: 400 / 1000 [ 40%]  (Warmup) 
#> Chain 1 Iteration: 500 / 1000 [ 50%]  (Warmup) 
#> Chain 1 Iteration: 501 / 1000 [ 50%]  (Sampling) 
#> Chain 1 Iteration: 600 / 1000 [ 60%]  (Sampling) 
#> Chain 1 Iteration: 700 / 1000 [ 70%]  (Sampling) 
#> Chain 1 Iteration: 800 / 1000 [ 80%]  (Sampling) 
#> Chain 1 Iteration: 900 / 1000 [ 90%]  (Sampling) 
#> Chain 1 Iteration: 1000 / 1000 [100%]  (Sampling) 
#> Chain 1 finished in 0.1 seconds.
#> Chain 2 Iteration:   1 / 1000 [  0%]  (Warmup) 
#> Chain 2 Iteration: 100 / 1000 [ 10%]  (Warmup) 
#> Chain 2 Iteration: 200 / 1000 [ 20%]  (Warmup) 
#> Chain 2 Iteration: 300 / 1000 [ 30%]  (Warmup) 
#> Chain 2 Iteration: 400 / 1000 [ 40%]  (Warmup) 
#> Chain 2 Iteration: 500 / 1000 [ 50%]  (Warmup) 
#> Chain 2 Iteration: 501 / 1000 [ 50%]  (Sampling) 
#> Chain 2 Iteration: 600 / 1000 [ 60%]  (Sampling) 
#> Chain 2 Iteration: 700 / 1000 [ 70%]  (Sampling) 
#> Chain 2 Iteration: 800 / 1000 [ 80%]  (Sampling) 
#> Chain 2 Iteration: 900 / 1000 [ 90%]  (Sampling) 
#> Chain 2 Iteration: 1000 / 1000 [100%]  (Sampling) 
#> Chain 2 finished in 0.1 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 0.1 seconds.
#> Total execution time: 0.4 seconds.
#> 
fit2 <- bsimm(
  ~Sex,
  mixture_data = sim$mixture_data,
  source_data = sim$source_data,
  tdf_data = sim$tdf_data,
  isotope_names = sim$isotope_names,
  source_means_sds = sim$source_means_sds,
  tdf_means_sds = sim$tdf_means_sds,
  conc_dep = sim$conc_dep,
  error_structure = sim$error_structure,
  source_col = sim$source_col,
  chains = 2,
  iter_warmup = 500,
  iter_sampling = 500
)
#> Running MCMC with 2 sequential chains...
#> 
#> Chain 1 Iteration:   1 / 1000 [  0%]  (Warmup) 
#> Chain 1 Iteration: 100 / 1000 [ 10%]  (Warmup) 
#> Chain 1 Iteration: 200 / 1000 [ 20%]  (Warmup) 
#> Chain 1 Iteration: 300 / 1000 [ 30%]  (Warmup) 
#> Chain 1 Iteration: 400 / 1000 [ 40%]  (Warmup) 
#> Chain 1 Iteration: 500 / 1000 [ 50%]  (Warmup) 
#> Chain 1 Iteration: 501 / 1000 [ 50%]  (Sampling) 
#> Chain 1 Iteration: 600 / 1000 [ 60%]  (Sampling) 
#> Chain 1 Iteration: 700 / 1000 [ 70%]  (Sampling) 
#> Chain 1 Iteration: 800 / 1000 [ 80%]  (Sampling) 
#> Chain 1 Iteration: 900 / 1000 [ 90%]  (Sampling) 
#> Chain 1 Iteration: 1000 / 1000 [100%]  (Sampling) 
#> Chain 1 finished in 0.3 seconds.
#> Chain 2 Iteration:   1 / 1000 [  0%]  (Warmup) 
#> Chain 2 Iteration: 100 / 1000 [ 10%]  (Warmup) 
#> Chain 2 Iteration: 200 / 1000 [ 20%]  (Warmup) 
#> Chain 2 Iteration: 300 / 1000 [ 30%]  (Warmup) 
#> Chain 2 Iteration: 400 / 1000 [ 40%]  (Warmup) 
#> Chain 2 Iteration: 500 / 1000 [ 50%]  (Warmup) 
#> Chain 2 Iteration: 501 / 1000 [ 50%]  (Sampling) 
#> Chain 2 Iteration: 600 / 1000 [ 60%]  (Sampling) 
#> Chain 2 Iteration: 700 / 1000 [ 70%]  (Sampling) 
#> Chain 2 Iteration: 800 / 1000 [ 80%]  (Sampling) 
#> Chain 2 Iteration: 900 / 1000 [ 90%]  (Sampling) 
#> Chain 2 Iteration: 1000 / 1000 [100%]  (Sampling) 
#> Chain 2 finished in 0.2 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 0.2 seconds.
#> Total execution time: 0.7 seconds.
#> 
loo::loo_compare(fit1, fit2)
#> Computing "loo" for model "fit1" (not cached; use `add_criterion()` to cache it
#> for reuse).
#> Computing "loo" for model "fit2" (not cached; use `add_criterion()` to cache it
#> for reuse).
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
#>  model elpd_diff se_diff p_worse diag_diff       diag_elpd
#>   fit2       0.0     0.0      NA           2 k_psis > 0.67
#>   fit1     -17.4     2.0    1.00   N < 100                
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
# }
```
