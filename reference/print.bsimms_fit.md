# Print a `bsimms` fit

Prints a one-screen overview of a fitted model: formula, source/isotope/
mixture counts and names, error structure, source/TDF data mode, whether
concentration dependence is enabled, group-level terms (if any), and the
fitting backend used. For posterior summaries, use
[`summary.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/summary.bsimms_fit.md).

## Usage

``` r
# S3 method for class 'bsimms_fit'
print(x, ...)
```

## Arguments

- x:

  A `bsimms_fit` object (as returned by
  [`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)).

- ...:

  Currently unused.

## Value

`x`, invisibly.

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
#> Chain 1 finished in 0.4 seconds.
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
#> Chain 2 finished in 0.3 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 0.3 seconds.
#> Total execution time: 0.9 seconds.
#> 
fit
#> Bayesian stable isotope mixing model (bsimms)
#>  formula:         ~1
#>  sources (K):     3 (Beaver, Deer, Hare)
#>  isotopes (J):    2 (d13C, d15N)
#>  mixtures (N):    10
#>  error structure: process_residual
#>  source data:     summary
#>  tdf data:        summary
#>  backend:         cmdstanr
# }
```
