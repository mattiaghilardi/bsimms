# Plot posterior source proportions

Plots posterior source proportions, as returned by
[`posterior_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/posterior_proportions.md)
or
[`fitted_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/fitted_proportions.md)
(`summary = FALSE`): a density or histogram of the posterior
distribution (`type = "density"`/`"histogram"`, one observation only),
or nested credible intervals across one or more observations
(`type = "interval"`, a forest/caterpillar plot coloured by source).
Requires `ggplot2`.

## Usage

``` r
plot_proportions(
  p_arr,
  type = c("density", "histogram", "interval"),
  probs = c(0.5, 0.95),
  robust = FALSE,
  point_size = 2,
  ...
)
```

## Arguments

- p_arr:

  A numeric `[n_draws, n_obs, K]` array of posterior proportion draws,
  with source names attached as the 3rd dimension's `dimnames`, as
  returned by
  [`posterior_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/posterior_proportions.md).

- type:

  Character; `"density"` (default), `"histogram"`, or `"interval"`.
  `"density"`/`"histogram"` require `p_arr` to have exactly one
  observation (row); use `"interval"` for more than one.

- probs:

  One or more credible-interval masses to display when
  `type = "interval"`, e.g. the default `c(0.5, 0.95)` draws both a 50%
  and a 95% interval. Ignored for `"density"`/`"histogram"`.

- robust:

  Logical; if `FALSE` (default) the point estimate (for
  `type = "interval"`) is the `mean`, if `TRUE` the `median`. Ignored
  for `"density"`/`"histogram"`.

- point_size:

  Size of the point marking the central estimate, for
  `type = "interval"` (default 2). Ignored otherwise.

- ...:

  Further arguments passed to the underlying `ggplot2` geom:
  [`ggplot2::geom_density()`](https://ggplot2.tidyverse.org/reference/geom_density.html),
  [`ggplot2::geom_histogram()`](https://ggplot2.tidyverse.org/reference/geom_histogram.html),
  or
  [`ggplot2::geom_linerange()`](https://ggplot2.tidyverse.org/reference/geom_linerange.html)
  for `type` `"density"`, `"histogram"`, or `"interval"` respectively.

## Value

A `ggplot` object.

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
#> Total execution time: 0.9 seconds.
#> 
p_arr <- posterior_proportions(fit)
plot_proportions(p_arr, type = "interval")

plot_proportions(p_arr[, 1, , drop = FALSE], type = "density")

# }
```
