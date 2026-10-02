# Approximate leave-one-out cross-validation

Computes PSIS-LOO (Vehtari, Gelman, and Gabry 2017) from the model's
`log_lik` generated quantity, which is the joint log density of each
mixture sample's full isotope profile (isotopes are summed/modelled
jointly per mixture sample, never held out individually, so the mixture
sample is the leave-one-out unit – see
[`make_stancode()`](https://mattiaghilardi.github.io/bsimms/reference/make_stancode.md)).

## Usage

``` r
# S3 method for class 'bsimms_fit'
loo(x, cores = getOption("mc.cores", 1), ...)
```

## Arguments

- x:

  A `bsimms_fit` object (as returned by
  [`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)).

- cores:

  Number of cores used for
  [`loo::relative_eff()`](https://mc-stan.org/loo/reference/relative_eff.html)
  and [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html)
  (default `getOption("mc.cores", 1)`).

- ...:

  Further arguments passed on to
  [`loo::loo.array()`](https://mc-stan.org/loo/reference/loo.html).

## Value

A `loo` object, as returned by
[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html).

## References

Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*, 27(5), 1413-1432.
[doi:10.1007/s11222-016-9696-4](https://doi.org/10.1007/s11222-016-9696-4)

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
#> Mean chain execution time: 0.4 seconds.
#> Total execution time: 0.9 seconds.
#> 
loo::loo(fit)
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
#> 
#> Computed from 1000 by 10 log-likelihood matrix.
#> 
#>          Estimate  SE
#> elpd_loo    -30.9 2.4
#> p_loo         3.1 0.7
#> looic        61.7 4.8
#> ------
#> MCSE of elpd_loo is NA.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.3, 0.6]).
#> 
#> Pareto k diagnostic values:
#>                           Count Pct.    Min. ESS
#> (-Inf, 0.67]   (good)     9     90.0%   71      
#>    (0.67, 1]   (bad)      1     10.0%   <NA>    
#>     (1, Inf)   (very bad) 0      0.0%   <NA>    
#> See help('pareto-k-diagnostic') for details.
# }
```
