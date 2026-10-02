# Cache LOO/WAIC/R-squared criteria on a `bsimms` fit

Computes one or more model-evaluation criteria and stores them in
`x$criteria` (a named list), returning the augmented fit.
[`loo_compare.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/loo_compare.bsimms_fit.md)
reuses a cached `"loo"`/`"waic"` criterion automatically instead of
recomputing it; `"bayes_R2"` is cached purely for reuse/comparison.

## Usage

``` r
add_criterion(x, ...)

# S3 method for class 'bsimms_fit'
add_criterion(x, criterion = "loo", ...)
```

## Arguments

- x:

  A `bsimms_fit` object (as returned by
  [`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)).

- ...:

  Further arguments passed on to
  [`loo.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/loo.bsimms_fit.md)
  (if `"loo"` is requested),
  [`waic.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/waic.bsimms_fit.md)
  (`"waic"`), or
  [`bayes_R2.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/bayes_R2.bsimms_fit.md)
  (`"bayes_R2"`, always cached with `summary = FALSE`, i.e. the raw
  posterior draws).

- criterion:

  Character vector of criteria to compute and cache: `"loo"` (default),
  `"waic"`, `"bayes_R2"`, or any combination.

## Value

`x`, with `x$criteria` updated to include the newly computed
criterion/criteria.

## Examples

``` r
# \donttest{
sim <- simulate_bsimms_data(
  ~1,
  n_mixture_obs = 10,
  source_names = c("Beaver", "Deer", "Hare"),
  isotope_names = c("d13C", "d15N"),
  source_means_sds = TRUE,
  seed = 1
)
fit <- bsimm(
  sim$formula,
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
#> Chain 1 finished in 0.5 seconds.
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
#> Chain 2 finished in 0.5 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 0.5 seconds.
#> Total execution time: 1.2 seconds.
#> 
fit <- add_criterion(fit, "loo")
fit$criteria$loo
#> 
#> Computed from 1000 by 10 log-likelihood matrix.
#> 
#>          Estimate  SE
#> elpd_loo    -30.7 2.2
#> p_loo         2.8 0.6
#> looic        61.4 4.4
#> ------
#> MCSE of elpd_loo is 0.1.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.3, 0.9]).
#> 
#> All Pareto k estimates are good (k < 0.67).
#> See help('pareto-k-diagnostic') for details.
# }
```
