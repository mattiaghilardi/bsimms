# Trace and density plots of model parameters

Plots the posterior density and MCMC trace of the model's underlying
Stan parameters (population-average source proportions `p_global`,
fixed-effect coefficients, group-level standard deviations, and error
term(s)), via
[`bayesplot::mcmc_combo()`](https://mc-stan.org/bayesplot/reference/MCMC-combos.html).
Use
[`plot_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/plot_proportions.md)
for summaries of the source proportions themselves, or
[`conditional_effects()`](https://mattiaghilardi.github.io/bsimms/reference/conditional_effects.md)
to see how they vary with a covariate.

## Usage

``` r
# S3 method for class 'bsimms_fit'
plot(
  x,
  variable = NULL,
  combo = c("dens", "trace"),
  nvariables = 5,
  plot = TRUE,
  ask = TRUE,
  newpage = TRUE,
  ...
)
```

## Arguments

- x:

  A `bsimms_fit` object (as returned by
  [`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)).

- variable:

  Optional character vector of parameter (base) names to plot. `NULL`
  (default) plots `p_global`, the fixed effects (if any), group-level
  standard deviations (if any), and the error term(s).

- combo:

  Character vector of two `bayesplot` `mcmc_*` plot types to combine
  (default `c("dens", "trace")`); see
  [`bayesplot::mcmc_combo()`](https://mc-stan.org/bayesplot/reference/MCMC-combos.html).

- nvariables:

  Maximum number of parameters shown per plot; models with more are
  split across multiple plots (default 5).

- plot:

  Logical; display each plot as a side effect (default `TRUE`)? If
  `FALSE`, the plots are only built and returned.

- ask:

  Logical; prompt the user before displaying each new page after the
  first (default `TRUE`)? Only relevant if `plot = TRUE` and there is
  more than one page.

- newpage:

  Logical; start the first page on a new graphics page (default `TRUE`)?
  Every page after the first always starts on a new page regardless.
  Only relevant if `plot = TRUE`.

- ...:

  Further arguments passed to
  [`bayesplot::mcmc_combo()`](https://mc-stan.org/bayesplot/reference/MCMC-combos.html).

## Value

A list of plot objects (one per page), invisibly; each is also displayed
as a side effect if `plot = TRUE`.

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
plot(fit)

# }
```
