# Simulate data for a `bsimms` model

Generates `mixture_data`, `source_data`, and `tdf_data` for an arbitrary
`bsimms` model configuration (any number of sources/isotopes, raw or
summarised source/TDF data, any error structure, with or without
concentration dependence, and any `lme4`-style fixed/random-effects
`formula` of mixture-level covariates), together with the known
generating ("true") parameter values. Useful for testing
[`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)
across configurations, and for checking parameter recovery.

## Usage

``` r
simulate_bsimms_data(
  formula = ~1,
  n_mixture_obs,
  n_sources = 3,
  source_names = NULL,
  n_isotopes = 2,
  isotope_names = NULL,
  n_levels = list(),
  n_groups = list(),
  balanced = TRUE,
  source_means_sds = FALSE,
  n_source_obs = 10,
  tdf_means_sds = TRUE,
  n_tdf_obs = 10,
  conc_dep = FALSE,
  error_structure = c("process_residual", "process_only", "residual_only"),
  p_global = NULL,
  sigma = NULL,
  resid_prop = NULL,
  source_col = "Source",
  seed = NULL
)
```

## Arguments

- formula:

  A one-sided `lme4`-style formula, as in
  [`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md).
  Random slopes (e.g. `(x | Group)`) are not supported.

- n_mixture_obs:

  Number of mixture observations to simulate.

- n_sources:

  Number of sources to simulate. Ignored if `source_names` is supplied.

- source_names:

  Optional character vector of source names. If `NULL` (default),
  `n_sources` sources named `"source1"`, `"source2"`, ... are used.

- n_isotopes:

  Number of isotopes to simulate. Ignored if `isotope_names` is
  supplied.

- isotope_names:

  Optional character vector of isotope names. If `NULL` (default),
  `n_isotopes` isotopes named `"isotope1"`, `"isotope2"`, ... are used.

- n_levels:

  Named list, one entry per fixed-effect factor in `formula`, giving
  that factor's number of levels (an integer \>= 2). Any variable in
  `formula` not listed here is generated as a continuous covariate
  instead.

- n_groups:

  Named list, one entry per random-effect grouping factor in `formula`
  (i.e. every variable referenced to the right of `|`), giving that
  factor's total number of groups (an integer \>= 2). Required for every
  grouping factor in `formula`. For a nested inner factor (e.g. `Site`
  in `(1 | Region/Site)`), this total is split across the outer factor's
  levels the same way `n_mixture_obs` is split across a factor's levels
  elsewhere (see `balanced`), so each outer level ends up with its own,
  possibly unequal, number of inner levels.

- balanced:

  Logical; split `n_mixture_obs` as evenly as possible across each
  factor's/grouping factor's levels (`TRUE`, default; sizes differ by at
  most one observation when `n_mixture_obs` is not a multiple of the
  number of levels), or assign them at random (`FALSE`)? Either way
  every level/group gets at least one observation.

- source_means_sds, tdf_means_sds:

  Logical; generate `source_data`/ `tdf_data` as means/SDs (`TRUE`) or
  raw replicate samples (`FALSE`)? Default `FALSE` for
  `source_means_sds`, `TRUE` for `tdf_means_sds`, matching
  [`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)'s
  own defaults.

- n_source_obs:

  Number of raw replicate observations to simulate per source: either a
  single integer, recycled to every source (default 10), or a numeric
  vector of length `n_sources` giving each source its own count
  (allowing unbalanced source data). Ignored if
  `source_means_sds = TRUE`.

- n_tdf_obs:

  As `n_source_obs`, for `tdf_data`. Ignored if `tdf_means_sds = TRUE`.

- conc_dep:

  Logical; also generate a `<isotope>_conc` column (in `(0, 1]`, summing
  to at most 1 across isotopes for a given source) in `source_data` for
  every isotope, enabling concentration dependence (default `FALSE`)?

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

- p_global:

  Optional numeric vector of length `n_sources`, the true baseline
  source proportions (must sum to 1). If `NULL` (default), randomly
  generated.

- sigma:

  Optional numeric vector of length `n_isotopes`, the true residual SD
  for each isotope. Used only when `error_structure = "residual_only"`
  (ignored otherwise). If `NULL` (default), randomly generated.

- resid_prop:

  Optional numeric vector of length `n_isotopes`, the true
  residual-error factor scaling process variance for each isotope. Used
  only when `error_structure = "process_residual"` (ignored otherwise).
  If `NULL` (default), randomly generated.

- source_col:

  Name of the source-identifier column shared by `source_data` and
  `tdf_data`. Default `"Source"`.

- seed:

  Optional integer seed. The global RNG state is saved and restored, so
  calling this function does not affect the caller's own random draws.

## Value

A list with elements `formula`, `mixture_data`, `source_data`,
`tdf_data`, `isotope_names`, `source_means_sds`, `tdf_means_sds`,
`conc_dep`, `error_structure`, `source_col` (in the same order as
[`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)'s
arguments, ready to pass to it), and `truth`, the true generating values
that
[`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)
would later estimate from the data: `p_global`, `fixed` (matrix of
fixed-effect ILR-space coefficients, if `formula` has fixed-effect
terms), `random` (list of group-level SDs/realised effects, one element
per random-effect term, if any), `source_mean`/`source_sd` (if
`source_means_sds = FALSE`), `tdf_mean`/`tdf_sd` (if
`tdf_means_sds = FALSE`), and `sigma` or `resid_prop` (whichever applies
to `error_structure`).

## Details

Source proportions and mixture isotope values are simulated forward from
known truth: a baseline `p_global`, formula-driven fixed/random-effect
deviations from it in ILR space, source and TDF isotope means/SDs, and
observation-level noise appropriate to `error_structure`. Random-effect
terms in `formula` are currently limited to random intercepts
(`(1 | Group)`, including crossed/nested grouping factors); random
slopes are not yet supported.

Every covariate referenced in `formula` must be either a fixed-effect
factor named in `n_levels`, a random-effect grouping factor named in
`n_groups`, or otherwise is generated as an independent standard normal
(continuous) covariate. Source means are generated to be well-separated
and non-collinear in isotope space (sampled without replacement from a
pool of `2 * n_sources` equally spaced candidate values per isotope,
redrawn if the resulting configuration is collinear), rather than
sampled independently from a probability distribution, which cannot
guarantee isotopically distinct sources.

A warning is issued when `n_sources > n_isotopes + 1`: the mixing system
is then underdetermined in the classical sense (more sources than
isotopes can resolve without help from the prior), which `bsimms`
supports but which relies more heavily on the prior than the data. A
warning is also issued for a nested pair whose inner total equals its
outer total, since that leaves exactly one inner group per outer group –
relabelling it rather than adding real nested replication.

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
  ~ Sex + (1 | Region),
  n_mixture_obs = 60,
  n_levels = list(Sex = 2),
  n_groups = list(Region = 3)
)
str(sim, max.level = 1)
#> List of 11
#>  $ formula         :Class 'formula'  language ~Sex + (1 | Region)
#>   .. ..- attr(*, ".Environment")=<environment: 0x5636863ae158> 
#>  $ mixture_data    :'data.frame':    60 obs. of  4 variables:
#>  $ source_data     :'data.frame':    30 obs. of  3 variables:
#>  $ tdf_data        :'data.frame':    3 obs. of  5 variables:
#>  $ isotope_names   : chr [1:2] "isotope1" "isotope2"
#>  $ source_means_sds: logi FALSE
#>  $ tdf_means_sds   : logi TRUE
#>  $ conc_dep        : logi FALSE
#>  $ error_structure : chr "process_residual"
#>  $ source_col      : chr "Source"
#>  $ truth           :List of 6
```
