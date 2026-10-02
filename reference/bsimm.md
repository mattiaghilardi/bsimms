# Fit a Bayesian stable isotope mixing model

The main entry point. Builds a bespoke Stan program for the requested
model (see
[`make_stancode()`](https://mattiaghilardi.github.io/bsimms/reference/make_stancode.md)),
compiles it, and draws posterior samples via MCMC (`cmdstanr` preferred,
`rstan` as a fallback).

## Usage

``` r
bsimm(
  formula,
  mixture_data,
  source_data,
  tdf_data,
  isotope_names,
  source_means_sds = FALSE,
  tdf_means_sds = TRUE,
  conc_dep = FALSE,
  error_structure = c("process_residual", "process_only", "residual_only"),
  prior = NULL,
  source_col = "Source",
  backend = c("auto", "cmdstanr", "rstan"),
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = NULL,
  cores = getOption("mc.cores", 1),
  refresh = NULL,
  ...
)
```

## Arguments

- formula:

  A one-sided `lme4`-style formula for the source proportions, e.g.
  `~ Sex + Season + (1 + Season | Region)`. All variables referenced
  must be columns of `mixture_data`. The population-level baseline is
  not part of `formula`: it is `p_global`, a `simplex[K]` parameter with
  a `Dirichlet(alpha)` prior (as in `MixSIAR`), and `formula`'s terms
  are covariate deviations from it, added in ILR space (see
  [`bsimms_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_prior.md)'s
  `"p_global"` class to set `alpha`). A formula with no covariates (e.g.
  `~ 1`) is a valid Dirichlet-baseline-only model.

- mixture_data:

  Data frame of mixture observations (e.g. consumers, in a diet-mixing
  study): one row per individual (or per sample), with a column for each
  entry of `isotope_names`, plus any covariate/grouping columns used in
  `formula`.

- source_data:

  Source isotope data. If `source_means_sds = FALSE` (default), long
  format with one row per raw sample and columns `Source`,
  `isotope_names`. If `source_means_sds = TRUE`, one row per source with
  columns `Source`, `<isotope>_mean`, `<isotope>_sd` for every isotope.

- tdf_data:

  Trophic discrimination factor (diet-tissue discrimination) data, same
  layout convention as `source_data`, controlled by `tdf_means_sds`
  (default `TRUE`, since TDFs are usually taken from the literature as
  means/SDs).

- isotope_names:

  Character vector of isotope column names shared by `mixture_data`,
  `source_data` and `tdf_data`, e.g. `c("d13C", "d15N")`.

- source_means_sds:

  Logical; is `source_data` supplied as means/SDs (`TRUE`) or raw
  replicate samples (`FALSE`, default)? When raw, source means and SDs
  are estimated as part of the model, and their uncertainty is fully
  propagated into the posterior.

- tdf_means_sds:

  Logical; is `tdf_data` supplied as means/SDs (`TRUE`, default) or raw
  replicate samples (`FALSE`)? When raw, TDF means and SDs are estimated
  as part of the model and their uncertainty is fully propagated into
  the posterior. Unlike raw source data, though, no cross-isotope
  correlation is estimated for raw TDF data (there is no `"tdf_cor"`
  class in
  [`bsimms_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_prior.md)):
  the mixture likelihood only takes the multivariate
  (`multi_normal_cholesky`) form when *source* data is raw and there are
  2+ isotopes. This asymmetry is intentional (raw TDF data with enough
  replicates per source to usefully estimate a correlation is a rare
  combination in practice), not an oversight, but could be added
  symmetrically to `source_cor` if needed.

- conc_dep:

  Logical; enable elemental concentration dependence (default `FALSE`)?
  If `TRUE`, `source_data` must have a `<isotope>_conc` column for every
  isotope in `isotope_names`, giving that source's proportional
  elemental concentration (in `(0, 1]`) for that isotope's element – one
  value per raw sample (averaged per source) if
  `source_means_sds = FALSE`, or one value per source if
  `source_means_sds = TRUE`. Since these are proportions of the same
  source's total mass, they cannot sum to more than 1 across isotopes.

- error_structure:

  One of `"process_residual"` (default: source/TDF variance propagated
  into the mixture and scaled by an estimated multiplicative
  residual-error factor, i.e. the "Residual \* Process" error structure
  of Stock & Semmens 2016), `"process_only"` (only source/TDF variance
  propagated, no separate residual term, as in the original MixSIR
  model, Moore & Semmens 2008), or `"residual_only"` (source/TDF
  variance is not propagated into the mixture at all; all unexplained
  variance instead goes into an isotope-specific residual error term, as
  in the original SIAR model, Parnell et al. 2010).

- prior:

  Optional `bsimms_prior` object (see
  [`bsimms_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_prior.md))
  with one or more rows overriding the default priors. Unspecified
  parameters keep their (weakly informative, partly data-scaled)
  defaults; see
  [`bsimms_get_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_get_prior.md).

- source_col:

  Name of the source-identifier column shared by `source_data` and
  `tdf_data`. Default `"Source"`.

- backend:

  `"cmdstanr"` or `"rstan"`. Default: `cmdstanr` if installed, otherwise
  `rstan`.

- chains:

  Positive integer, number of MCMC chains (default 4).

- iter_warmup:

  Positive integer, warmup iterations per chain (default 1000).

- iter_sampling:

  Positive integer, post-warmup sampling iterations per chain (default
  1000).

- seed:

  Optional integer seed.

- cores:

  Positive integer, number of cores for parallel chains (default
  `getOption("mc.cores", 1)`).

- refresh:

  How often to print sampler progress (in iterations); default lets the
  backend choose.

- ...:

  Further arguments passed on to `cmdstanr`'s `$sample()` or
  [`rstan::sampling()`](https://mc-stan.org/rstan/reference/stanmodel-method-sampling.html).

## Value

An object of class `bsimms_fit`: a list with elements `fit` (the raw
backend fit object), `backend`, `spec` (internal model specification),
`stancode`, `standata`, `prior`, and `call`.

## References

Moore, J.W., & Semmens, B.X. (2008). Incorporating uncertainty and prior
information into stable isotope mixing models. *Ecology Letters*, 11(5),
470-480.
[doi:10.1111/j.1461-0248.2008.01163.x](https://doi.org/10.1111/j.1461-0248.2008.01163.x)

Parnell, A.C., Inger, R., Bearhop, S., & Jackson, A.L. (2010). Source
partitioning using stable isotopes: coping with too much variation.
*PLoS ONE*, 5(3), e9672.
[doi:10.1371/journal.pone.0009672](https://doi.org/10.1371/journal.pone.0009672)

Stock, B.C., & Semmens, B.X. (2016). Unifying error structures in
commonly used biotracer mixing models. *Ecology*, 97(10), 2562-2569.
[doi:10.1002/ecy.1517](https://doi.org/10.1002/ecy.1517)

## Examples

``` r
# \donttest{
sim <- simulate_bsimms_data(
  ~1,
  n_mixture_obs = 10,
  source_names = c("Beaver", "Deer", "Hare"),
  isotope_names = c("d13C", "d15N"),
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
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: cholesky_decompose: A is not symmetric. A[1,2] = -nan, but A[2,1] = -nan (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 102, column 4 to column 83)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: cholesky_decompose: A is not symmetric. A[1,2] = -nan, but A[2,1] = -nan (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 102, column 4 to column 83)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: cholesky_decompose: A is not symmetric. A[1,2] = -nan, but A[2,1] = -nan (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 102, column 4 to column 83)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: cholesky_decompose: A is not symmetric. A[1,2] = -nan, but A[2,1] = -nan (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 102, column 4 to column 83)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: cholesky_decompose: A is not symmetric. A[1,2] = -inf, but A[2,1] = -inf (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 102, column 4 to column 83)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lkj_corr_cholesky_lpdf: Random variable[2] is 0, but must be positive! (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 129, column 2 to column 41)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lkj_corr_cholesky_lpdf: Random variable[2] is 0, but must be positive! (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 130, column 2 to column 41)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lkj_corr_cholesky_lpdf: Random variable[2] is 0, but must be positive! (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 130, column 2 to column 41)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: cholesky_decompose: A is not symmetric. A[1,2] = -nan, but A[2,1] = -nan (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 102, column 4 to column 83)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: cholesky_decompose: A is not symmetric. A[1,2] = -nan, but A[2,1] = -nan (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 102, column 4 to column 83)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
#> Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 1 Exception: lkj_corr_cholesky_lpdf: Random variable[2] is 0, but must be positive! (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 129, column 2 to column 41)
#> Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 1 
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
#> Chain 1 finished in 1.7 seconds.
#> Chain 2 Iteration:   1 / 1000 [  0%]  (Warmup) 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: cholesky_decompose: A is not symmetric. A[1,2] = -nan, but A[2,1] = -nan (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 102, column 4 to column 83)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: cholesky_decompose: A is not symmetric. A[1,2] = -nan, but A[2,1] = -nan (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 102, column 4 to column 83)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: cholesky_decompose: A is not symmetric. A[1,2] = -nan, but A[2,1] = -nan (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 102, column 4 to column 83)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: cholesky_decompose: A is not symmetric. A[1,2] = inf, but A[2,1] = inf (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 102, column 4 to column 83)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lkj_corr_cholesky_lpdf: Random variable[2] is 0, but must be positive! (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 129, column 2 to column 41)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lkj_corr_cholesky_lpdf: Random variable[2] is 0, but must be positive! (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 129, column 2 to column 41)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: cholesky_decompose: A is not symmetric. A[1,2] = -nan, but A[2,1] = -nan (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 102, column 4 to column 83)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: cholesky_decompose: A is not symmetric. A[1,2] = -nan, but A[2,1] = -nan (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 102, column 4 to column 83)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
#> Chain 2 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
#> Chain 2 Exception: lkj_corr_cholesky_lpdf: Random variable[2] is 0, but must be positive! (in '/tmp/RtmpfMAX2p/model-1b8c77830c2.stan', line 129, column 2 to column 41)
#> Chain 2 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
#> Chain 2 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
#> Chain 2 
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
#> Chain 2 finished in 1.9 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 1.8 seconds.
#> Total execution time: 3.7 seconds.
#> 
# }
```
