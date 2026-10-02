# Default priors for a `bsimms` model

Builds the model specification (design matrices, source/TDF data layout)
without generating or fitting the Stan model, and returns the table of
default priors that would be used. Edit rows of the result and pass to
the `prior` argument of
[`make_stancode()`](https://mattiaghilardi.github.io/bsimms/reference/make_stancode.md),
[`make_standata()`](https://mattiaghilardi.github.io/bsimms/reference/make_standata.md)
or
[`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)
to override.

## Usage

``` r
bsimms_get_prior(
  formula,
  mixture_data,
  source_data,
  tdf_data,
  isotope_names,
  source_means_sds = FALSE,
  tdf_means_sds = TRUE,
  conc_dep = FALSE,
  error_structure = c("process_residual", "process_only", "residual_only"),
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

- source_col:

  Name of the source-identifier column shared by `source_data` and
  `tdf_data`. Default `"Source"`.

## Value

A `bsimms_prior` data frame.

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
[doi:10.1002/ecy.1517](https://doi.org/10.1002/ecy.1517) Appendix S2
simulation-tests the default `"resid_prop"` prior (`uniform(0, 20)`)
against three priors centred/moded at 1 (chi-square(3), gamma(2,2),
log-normal(0,1)) and finds it better-calibrated at high consumer
variance.

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
bsimms_get_prior(
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
#>                      prior       class coef resp  group
#>               normal(0, 1)           b                 
#>                          1    p_global           Beaver
#>                          1    p_global             Deer
#>                          1    p_global             Hare
#>         student_t(3, 0, 1)          sd                 
#>       lkj_corr_cholesky(1)         cor                 
#>             uniform(0, 20)  resid_prop                 
#>             uniform(0, 20)  resid_prop      d13C       
#>             uniform(0, 20)  resid_prop      d15N       
#>   normal(2.65907, 7.01714) source_mean      d13C Beaver
#>  student_t(3, 0, 0.701714)   source_sd      d13C Beaver
#>  normal(-7.62958, 7.22181) source_mean      d13C   Deer
#>  student_t(3, 0, 0.722181)   source_sd      d13C   Deer
#>   normal(12.4633, 7.77402) source_mean      d13C   Hare
#>  student_t(3, 0, 0.777402)   source_sd      d13C   Hare
#>  normal(-12.3577, 7.73755) source_mean      d15N Beaver
#>  student_t(3, 0, 0.773755)   source_sd      d15N Beaver
#>   normal(12.5269, 6.85126) source_mean      d15N   Deer
#>  student_t(3, 0, 0.685126)   source_sd      d15N   Deer
#>   normal(1.46555, 3.66305) source_mean      d15N   Hare
#>  student_t(3, 0, 0.366305)   source_sd      d15N   Hare
#>       lkj_corr_cholesky(1)  source_cor                 
#>       lkj_corr_cholesky(1)  source_cor           Beaver
#>       lkj_corr_cholesky(1)  source_cor             Deer
#>       lkj_corr_cholesky(1)  source_cor             Hare
```
