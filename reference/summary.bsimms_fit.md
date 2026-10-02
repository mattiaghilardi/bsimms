# Summarise a `bsimms` fit

Reports posterior summaries for the global (population-average) source
proportions (`p_global`, the `Dirichlet`-distributed baseline shared by
all mixture samples, as in `MixSIAR`), the ILR-scale fixed-effect
covariate slopes (deviations from that baseline, if any), group-level
standard deviations (if any), and the error term(s).

## Usage

``` r
# S3 method for class 'bsimms_fit'
summary(object, robust = FALSE, probs = c(0.025, 0.975), ...)
```

## Arguments

- object:

  A `bsimms_fit` object (as returned by
  [`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)).

- robust:

  Logical; if `FALSE` (default) summarise the central tendency/spread
  with `mean`/`sd`, if `TRUE` use `median`/`mad` instead.

- probs:

  Quantiles to report. Default 2.5%/97.5%.

- ...:

  Currently unused.

## Value

An object of class `summary.bsimms_fit`.

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
summary(fit)
#> Bayesian stable isotope mixing model
#>  formula: ~1
#> 
#> Population-average source proportions:
#>  source  mean    sd  q2.5 q97.5  rhat ess_bulk ess_tail
#>  Beaver 0.208 0.106 0.036 0.442 1.014      320      215
#>    Deer 0.393 0.065 0.258 0.503 1.011      365      397
#>    Hare 0.400 0.045 0.296 0.474 1.013      321      333
#> 
#> Error term(s):
#>  isotope  mean    sd  q2.5 q97.5  rhat ess_bulk ess_tail
#>     d13C 3.735 2.397 1.177 9.637 1.002      511      272
#>     d15N 1.801 1.183 0.587 4.540 1.000      481      459
# }
```
