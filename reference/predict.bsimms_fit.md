# Posterior predictive draws of mixture isotope values (summarised)

Returns posterior summaries of the model's posterior predictive draws
(i.e.
[`fitted.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/fitted.bsimms_fit.md)
plus observation error), via
[`posterior_predict.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/posterior_predict.bsimms_fit.md).

## Usage

``` r
# S3 method for class 'bsimms_fit'
predict(
  object,
  newdata = NULL,
  resp = NULL,
  re_formula = NULL,
  allow_new_levels = FALSE,
  sample_new_levels = c("uncertainty", "gaussian"),
  ndraws = NULL,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975),
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

- summary:

  Logical; return posterior summaries (default) or raw draws (`FALSE`,
  equivalent to
  [`posterior_predict.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/posterior_predict.bsimms_fit.md)).

- robust:

  Logical; if `FALSE` (default) summarise the central tendency/spread
  with `mean`/`sd`, if `TRUE` use `median`/`mad` instead.

- probs:

  Quantiles to report when `summary = TRUE`.

- ...:

  Currently unused.

## Value

If `summary = TRUE`, a long-format data frame (`row`, `isotope`, and one
column per summary measure, named as in
[`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html));
if `FALSE`, an array `[n_draws, n_obs, length(resp)]`.

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
#> Chain 2 finished in 0.2 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 0.2 seconds.
#> Total execution time: 0.6 seconds.
#> 
predict(fit)
#>    row isotope       mean       sd       q2.5    q97.5
#> 1    1    d13C  3.6468671 1.510198  0.6091445 6.502106
#> 2    1    d15N -0.3740708 1.118937 -2.5655796 1.800165
#> 3    2    d13C  3.5111371 1.445808  0.5620106 6.274583
#> 4    2    d15N -0.4668943 1.055497 -2.5823317 1.692926
#> 5    3    d13C  3.6199466 1.435685  0.6575557 6.402304
#> 6    3    d15N -0.4411145 1.094167 -2.7021961 1.690683
#> 7    4    d13C  3.5985799 1.434113  0.6336140 6.293999
#> 8    4    d15N -0.4294275 1.114727 -2.6318901 1.974058
#> 9    5    d13C  3.5769260 1.580797  0.4369639 6.651124
#> 10   5    d15N -0.4311221 1.077828 -2.5895067 1.648285
#> 11   6    d13C  3.6459406 1.451426  0.9244438 6.557090
#> 12   6    d15N -0.4764442 1.093526 -2.6313774 1.790522
#> 13   7    d13C  3.5809305 1.427078  0.6963427 6.358213
#> 14   7    d15N -0.4546859 1.072438 -2.5491880 1.710965
#> 15   8    d13C  3.6055988 1.461097  0.6060846 6.263207
#> 16   8    d15N -0.3704034 1.056319 -2.6448879 1.958759
#> 17   9    d13C  3.5705164 1.506772  0.5243990 6.435281
#> 18   9    d15N -0.4122829 1.062361 -2.6074980 1.620477
#> 19  10    d13C  3.5738800 1.442315  0.6037786 6.400365
#> 20  10    d15N -0.4063002 1.049020 -2.5045544 1.784974
# }
```
