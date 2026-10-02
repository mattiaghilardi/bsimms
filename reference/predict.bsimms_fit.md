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
#> Total execution time: 1.1 seconds.
#> 
predict(fit)
#>    row isotope       mean        sd       q2.5    q97.5
#> 1    1    d13C  3.6278100 1.5868531  0.4892001 6.612824
#> 2    1    d15N -0.4319178 1.0336033 -2.5030207 1.590887
#> 3    2    d13C  3.5781958 1.4845521  0.6451600 6.672481
#> 4    2    d15N -0.4507872 1.0507375 -2.6502427 1.512005
#> 5    3    d13C  3.6222577 1.5534812  0.2729317 6.798991
#> 6    3    d15N -0.4742205 1.0819579 -2.7111523 1.663810
#> 7    4    d13C  3.6091725 1.4827827  0.6351185 6.465026
#> 8    4    d15N -0.4466397 0.9937387 -2.4683207 1.399873
#> 9    5    d13C  3.6339085 1.5420587  0.4752791 6.660485
#> 10   5    d15N -0.4294893 1.0484032 -2.5084212 1.668854
#> 11   6    d13C  3.5948796 1.5761043  0.5725467 6.679204
#> 12   6    d15N -0.4190906 1.0659680 -2.5812912 1.638024
#> 13   7    d13C  3.5116354 1.5038740  0.2781620 6.502530
#> 14   7    d15N -0.4812810 1.0430870 -2.6255353 1.642433
#> 15   8    d13C  3.5684010 1.5282437  0.4602685 6.396085
#> 16   8    d15N -0.4208288 1.1005208 -2.5885984 1.734753
#> 17   9    d13C  3.6403469 1.5853255  0.2361198 6.854503
#> 18   9    d15N -0.4573589 0.9804324 -2.2775902 1.508166
#> 19  10    d13C  3.5963131 1.5234963  0.2380169 6.681450
#> 20  10    d15N -0.4525511 1.0503825 -2.5795652 1.740258
# }
```
