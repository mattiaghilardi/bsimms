# Generate the Stan data list for a `bsimms` model

Builds the named list of data passed to Stan, matching exactly the
`data` block produced by
[`make_stancode()`](https://mattiaghilardi.github.io/bsimms/reference/make_stancode.md)
for the same arguments.

## Usage

``` r
make_standata(
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
  source_col = "Source"
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

## Value

A named list suitable for `cmdstanr::sample(data = ...)` or
`rstan::sampling(data = ...)`. Class `bsimms_standata`.

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
sim <- simulate_bsimms_data(
  ~ 1 + (1 | Region),
  n_mixture_obs = 20,
  source_names = c("Beaver", "Deer", "Hare"),
  isotope_names = c("d13C", "d15N"),
  n_groups = list(Region = 2),
  seed = 1
)
sdata <- make_standata(
  sim$formula,
  mixture_data = sim$mixture_data,
  source_data = sim$source_data,
  tdf_data = sim$tdf_data,
  isotope_names = sim$isotope_names,
  source_means_sds = sim$source_means_sds,
  tdf_means_sds = sim$tdf_means_sds,
  conc_dep = sim$conc_dep,
  error_structure = sim$error_structure,
  source_col = sim$source_col
)
str(sdata, max.level = 1)
#> List of 16
#>  $ N              : int 20
#>  $ J              : int 2
#>  $ K              : int 3
#>  $ D              : int 2
#>  $ V              : num [1:3, 1:2] 0.707 -0.707 0 0.408 0.408 ...
#>  $ y              : num [1:20, 1:2] 3.76 3.05 2.44 2.22 2.5 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ alpha_dirichlet: num [1:3] 1 1 1
#>  $ N_re_Region    : int 2
#>  $ M_re_Region    : int 1
#>  $ Z_re_Region    : num [1:20, 1] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   ..- attr(*, "assign")= int 0
#>  $ grp_re_Region  : int [1:20] 1 1 1 1 2 2 2 2 2 1 ...
#>  $ N_source_raw   : int 30
#>  $ source_idx     : int [1:30] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ source_raw     : num [1:30, 1:2] 1.91 3.43 3.08 2.34 2.18 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ tdf_mean_data  : num [1:3, 1:2] 1.95 1.89 1.91 1.87 2.15 ...
#>  $ tdf_sd_data    : num [1:3, 1:2] 0.636 0.439 0.434 0.491 0.657 ...
#>  - attr(*, "class")= chr [1:2] "bsimms_standata" "list"
```
