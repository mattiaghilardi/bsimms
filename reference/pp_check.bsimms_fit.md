# Posterior predictive checks

Compares observed mixture isotope values against posterior predictive
draws (`y_rep`) using
[`bayesplot::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html)'s
`ppc_*` plotting functions. Since these expect a single response vector,
one isotope must be selected.

## Usage

``` r
# S3 method for class 'bsimms_fit'
pp_check(object, resp = NULL, type = "dens_overlay", ndraws = NULL, ...)
```

## Arguments

- object:

  A `bsimms_fit` object (as returned by
  [`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)).

- resp:

  Character; the isotope to plot (one of `isotope_names`). Required if
  the model has more than one isotope; defaults to the only isotope
  otherwise.

- type:

  Character; the `bayesplot` PPC plot type, e.g. `"dens_overlay"`
  (default), `"hist"`, `"stat"`, `"scatter_avg"`, `"intervals"` — see
  [`bayesplot::available_ppc()`](https://mc-stan.org/bayesplot/reference/available_ppc.html)
  for the full list (passed without its `"ppc_"` prefix).

- ndraws:

  Optional integer; number of `y_rep` draws to (randomly) subsample for
  the plot. `NULL` (default) uses every draw for PPC types that
  aggregate/summarise across draws (e.g. `"stat"`, `"intervals"`,
  `"scatter_avg"`), where more draws only improve precision, or 10 draws
  for types that overlay one `y_rep` dataset per draw (e.g.
  `"dens_overlay"`, `"hist"`), where using every draw would overplot the
  figure.

- ...:

  Further arguments passed on to the underlying `ppc_*` function, e.g.
  `group` for grouped types.

## Value

A ggplot object, as returned by the underlying `ppc_*` function.

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
#> Chain 2 finished in 0.4 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 0.4 seconds.
#> Total execution time: 1.0 seconds.
#> 
bayesplot::pp_check(fit, resp = "d13C")
#> Using 10 posterior draws for ppc type "dens_overlay" by default.

# }
```
