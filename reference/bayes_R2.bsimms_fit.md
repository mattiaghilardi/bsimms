# Bayesian R-squared

Computes a Bayesian R-squared (Gelman, Goodrich, Gabry, and Vehtari
2019) for one or more isotopes. If `object$criteria$bayes_R2` was
already cached via
[`add_criterion()`](https://mattiaghilardi.github.io/bsimms/reference/add_criterion.md),
it is reused (subset to `resp`) instead of recomputed.

## Usage

``` r
# S3 method for class 'bsimms_fit'
bayes_R2(
  object,
  resp = NULL,
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

- resp:

  Optional character vector of isotope names (a subset of
  `isotope_names`) to compute R-squared for. `NULL` (default) computes
  it for every isotope.

- summary:

  Logical; return a summary data frame (default) or the raw posterior
  draws.

- robust:

  Logical; if `FALSE` (default) summarise the central tendency/spread
  with `mean`/`sd`, if `TRUE` use `median`/`mad` instead.

- probs:

  Quantiles to include in the summary (default 2.5%/97.5%).

- ...:

  Currently unused.

## Value

If `summary = TRUE`, a data frame with one row per isotope (`isotope`,
and one column per summary measure, named as in
[`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html):
`mean`, `sd`, `q2.5`, `q97.5`, or with `robust = TRUE`, `median`, `mad`,
`q2.5`, `q97.5`). If `FALSE`, an `n_draws x length(resp)` matrix of
R-squared draws (column names = isotope names).

## References

Gelman, A., Goodrich, B., Gabry, J., & Vehtari, A. (2019). R-squared for
Bayesian regression models. *The American Statistician*, 73(3), 307-309.
[doi:10.1080/00031305.2018.1549100](https://doi.org/10.1080/00031305.2018.1549100)
.
([Preprint](https://acris.aalto.fi/ws/portalfiles/portal/34206843/bayes_R2_v3.pdf))

## Examples

``` r
# \donttest{
sim <- simulate_bsimms_data(
  ~Sex,
  n_mixture_obs = 10,
  source_names = c("Beaver", "Deer", "Hare"),
  isotope_names = c("d13C", "d15N"),
  source_means_sds = TRUE,
  n_levels = list(Sex = 2),
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
#> Chain 1 finished in 0.2 seconds.
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
rstantools::bayes_R2(fit)
#>   isotope      mean         sd       q2.5     q97.5
#> 1    d13C 0.9510045 0.03618861 0.89811626 0.9638898
#> 2    d15N 0.6087731 0.17573367 0.08119061 0.7568052
# }
```
