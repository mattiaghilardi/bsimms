# Extract posterior draws from a `bsimms` fit

Backend-agnostic accessor returning a
[`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html)
regardless of whether the model was fit with `cmdstanr` or `rstan`.

## Usage

``` r
bsimms_draws(object, variable = NULL)
```

## Arguments

- object:

  A `bsimms_fit` object (as returned by
  [`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)).

- variable:

  Optional character vector of variable names (or `posterior`-style
  selectors) to extract; `NULL` extracts everything.

## Value

A
[`posterior::draws_array`](https://mc-stan.org/posterior/reference/draws_array.html).

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
#> Total execution time: 1.0 seconds.
#> 
bsimms_draws(fit, variable = "p_global")
#> # A draws_array: 500 iterations, 2 chains, and 3 variables
#> , , variable = p_global[1]
#> 
#>          chain
#> iteration    1    2
#>         1 0.19 0.31
#>         2 0.25 0.30
#>         3 0.22 0.30
#>         4 0.21 0.27
#>         5 0.19 0.33
#> 
#> , , variable = p_global[2]
#> 
#>          chain
#> iteration    1    2
#>         1 0.41 0.33
#>         2 0.36 0.36
#>         3 0.39 0.33
#>         4 0.38 0.38
#>         5 0.41 0.31
#> 
#> , , variable = p_global[3]
#> 
#>          chain
#> iteration    1    2
#>         1 0.40 0.36
#>         2 0.39 0.34
#>         3 0.39 0.37
#>         4 0.41 0.36
#>         5 0.40 0.35
#> 
#> # ... with 495 more iterations
# }
```
