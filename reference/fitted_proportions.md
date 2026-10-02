# Posterior source proportions (summarised)

Returns posterior summaries of source proportions `p`, via
[`posterior_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/posterior_proportions.md),
either for the mixture samples in the fitted data (`newdata = NULL`) or
for new mixture covariate combinations.

## Usage

``` r
fitted_proportions(
  object,
  newdata = NULL,
  re_formula = NULL,
  allow_new_levels = FALSE,
  sample_new_levels = c("uncertainty", "gaussian"),
  ndraws = NULL,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975)
)
```

## Arguments

- object:

  A `bsimms_fit` object (as returned by
  [`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)).

- newdata:

  Optional data frame of new mixture covariate values (same columns as
  used in `formula`).

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
  [`posterior_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/posterior_proportions.md)).

- robust:

  Logical; if `FALSE` (default) summarise the central tendency/spread
  with `mean`/`sd`, if `TRUE` use `median`/`mad` instead.

- probs:

  Quantiles to report when `summary = TRUE`.

## Value

If `summary = TRUE`, a long-format data frame with one row per
(observation, source): columns `row`, `source`, and one column per
summary measure, named as in
[`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html)
(`mean`, `sd`, `q2.5`, `q97.5`, or with `robust = TRUE`, `median`,
`mad`, `q2.5`, `q97.5`). If `summary = FALSE`, a numeric array
`[n_draws, n_obs, K]`.

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
fitted_proportions(fit)
#>    row source      mean         sd       q2.5     q97.5
#> 1    1 Beaver 0.2002756 0.10320212 0.02201473 0.4216104
#> 2    1   Deer 0.3960742 0.06266644 0.26609579 0.5051173
#> 3    1   Hare 0.4036503 0.04371577 0.31031700 0.4814764
#> 4    2 Beaver 0.2002756 0.10320212 0.02201473 0.4216104
#> 5    2   Deer 0.3960742 0.06266644 0.26609579 0.5051173
#> 6    2   Hare 0.4036503 0.04371577 0.31031700 0.4814764
#> 7    3 Beaver 0.2002756 0.10320212 0.02201473 0.4216104
#> 8    3   Deer 0.3960742 0.06266644 0.26609579 0.5051173
#> 9    3   Hare 0.4036503 0.04371577 0.31031700 0.4814764
#> 10   4 Beaver 0.2002756 0.10320212 0.02201473 0.4216104
#> 11   4   Deer 0.3960742 0.06266644 0.26609579 0.5051173
#> 12   4   Hare 0.4036503 0.04371577 0.31031700 0.4814764
#> 13   5 Beaver 0.2002756 0.10320212 0.02201473 0.4216104
#> 14   5   Deer 0.3960742 0.06266644 0.26609579 0.5051173
#> 15   5   Hare 0.4036503 0.04371577 0.31031700 0.4814764
#> 16   6 Beaver 0.2002756 0.10320212 0.02201473 0.4216104
#> 17   6   Deer 0.3960742 0.06266644 0.26609579 0.5051173
#> 18   6   Hare 0.4036503 0.04371577 0.31031700 0.4814764
#> 19   7 Beaver 0.2002756 0.10320212 0.02201473 0.4216104
#> 20   7   Deer 0.3960742 0.06266644 0.26609579 0.5051173
#> 21   7   Hare 0.4036503 0.04371577 0.31031700 0.4814764
#> 22   8 Beaver 0.2002756 0.10320212 0.02201473 0.4216104
#> 23   8   Deer 0.3960742 0.06266644 0.26609579 0.5051173
#> 24   8   Hare 0.4036503 0.04371577 0.31031700 0.4814764
#> 25   9 Beaver 0.2002756 0.10320212 0.02201473 0.4216104
#> 26   9   Deer 0.3960742 0.06266644 0.26609579 0.5051173
#> 27   9   Hare 0.4036503 0.04371577 0.31031700 0.4814764
#> 28  10 Beaver 0.2002756 0.10320212 0.02201473 0.4216104
#> 29  10   Deer 0.3960742 0.06266644 0.26609579 0.5051173
#> 30  10   Hare 0.4036503 0.04371577 0.31031700 0.4814764
# }
```
