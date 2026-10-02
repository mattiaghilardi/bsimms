# Reshape a draws array into long format

Reshapes a `[n_draws, n_obs, n_var]` array (as returned by
[`posterior_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/posterior_proportions.md),
[`posterior_epred.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/posterior_epred.bsimms_fit.md),
or
[`posterior_predict.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/posterior_predict.bsimms_fit.md)),
with variable names attached as the 3rd dimension's `dimnames`, into a
long-format data frame with one row per (draw, observation, variable)
triple, ready for custom plots or summaries (e.g. with
`ggplot2`/`ggdist`).

## Usage

``` r
draws_long(arr, var_col = "variable", value_col = "value")
```

## Arguments

- arr:

  A numeric `[n_draws, n_obs, n_var]` array, with variable names
  attached as the 3rd dimension's `dimnames`.

- var_col:

  Name to give the variable-identity column (default `"variable"`), e.g.
  `"source"` for a
  [`posterior_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/posterior_proportions.md)
  array or `"isotope"` for a
  [`posterior_epred.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/posterior_epred.bsimms_fit.md)
  array.

- value_col:

  Name to give the value column (default `"value"`).

## Value

A long-format data frame with columns `draw` (posterior draw index),
`row` (observation index, into `newdata` or the fitted mixture samples),
`<var_col>`, and `<value_col>`.

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
#> Chain 2 finished in 0.4 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 0.4 seconds.
#> Total execution time: 1.0 seconds.
#> 
p_arr <- posterior_proportions(fit)
long <- draws_long(p_arr, var_col = "source", value_col = "proportion")
head(long)
#>   draw row source proportion
#> 1    1   1 Beaver 0.20041774
#> 2    2   1 Beaver 0.28020488
#> 3    3   1 Beaver 0.06652705
#> 4    4   1 Beaver 0.20382201
#> 5    5   1 Beaver 0.16805357
#> 6    6   1 Beaver 0.18870313
# }
```
