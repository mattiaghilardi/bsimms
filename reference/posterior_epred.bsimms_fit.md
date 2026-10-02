# Expected mixture isotope values

Returns raw posterior draws of the model's expected mixture isotope
value `mu` (the mixing-model prediction *before* observation error is
added), either for the fitted mixture samples (`newdata = NULL`) or for
new mixture covariate combinations. Use
[`fitted.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/fitted.bsimms_fit.md)
for a summarised version instead,
[`posterior_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/posterior_proportions.md)
for the underlying source proportions, or
[`posterior_predict.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/posterior_predict.bsimms_fit.md)
for posterior predictive draws that include observation error.

## Usage

``` r
# S3 method for class 'bsimms_fit'
posterior_epred(
  object,
  newdata = NULL,
  resp = NULL,
  re_formula = NULL,
  allow_new_levels = FALSE,
  sample_new_levels = c("uncertainty", "gaussian"),
  ndraws = NULL,
  ...
)
```

## Arguments

- object:

  A `bsimms_fit` object (as returned by
  [`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)).

- newdata:

  Optional data frame of new mixture covariate values (same columns as
  used in `formula`).

- resp:

  Character vector of isotope name(s) (from `isotope_names`) to return.
  `NULL` (default) returns all isotopes.

- re_formula:

  Which group-level (random-effect) terms to condition on, whether
  predicting for the fitted mixture samples (`newdata = NULL`) or for
  new ones: `NULL` (default) conditions on every group-level term in the
  fitted model (with `newdata`, every term's grouping column(s) must
  then be supplied); `NA` or `~0` conditions on none of them
  (population-average for every term, regardless of what `newdata`
  contains); a reduced formula naming a subset, e.g. `~ (1 | Site)`,
  conditions only on the named term(s) (with `newdata`, only those
  term(s)' columns are needed). When `newdata` is `NULL`,
  `allow_new_levels`/`sample_new_levels` are ignored, since a fitted
  mixture sample's group membership is never new.

- allow_new_levels:

  Logical; only relevant with `newdata`. If `FALSE` (default), a
  grouping column value not seen when fitting (including `NA`, which is
  always treated as a new level) raises an error. If `TRUE`, a posterior
  draw is instead sampled for that new level via `sample_new_levels`;
  every row sharing the same new, non-`NA` level value gets the same
  sampled draws (as for a real level), while each `NA` row is sampled
  independently (it asserts no shared identity).

- sample_new_levels:

  How to sample a new level's group-level deviation when
  `allow_new_levels = TRUE`: `"uncertainty"` (default) draws, for each
  posterior draw, from a randomly chosen *existing* level's draw at that
  same iteration – most appropriate with many existing levels, where
  this empirically reflects the observed between-level spread;
  `"gaussian"` instead draws a fresh value for the new level from a
  normal distribution centred at zero, using that draw's own estimated
  group-level SD (and correlation, for multi-column terms) – more
  appropriate with few existing levels, where `"uncertainty"` could only
  ever resample one of a handful of specific levels rather than
  represent a genuinely new one.

- ndraws:

  Number of posterior draws to use, randomly subset from the full
  posterior. `NULL` (default) uses all draws.

- ...:

  Currently unused.

## Value

A numeric `[n_draws, n_obs, length(resp)]` array, with isotope names
attached as the 3rd dimension's `dimnames`.

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
#> Chain 2 finished in 0.5 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 0.4 seconds.
#> Total execution time: 1.1 seconds.
#> 
mu_arr <- rstantools::posterior_epred(fit)
dim(mu_arr)
#> [1] 1000   10    2
# }
```
