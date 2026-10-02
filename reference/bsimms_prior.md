# Specify a prior for one or more `bsimms` parameters

Build up a full prior specification by combining several calls with
[`c()`](https://rdrr.io/r/base/c.html), e.g.
`c(bsimms_prior("normal(0, 2)", class = "b"), bsimms_prior("normal(1, 1)", class = "b", coef = "SeasonWinter"))`.
More specific rows (a given `coef`/`resp`/`group`) take precedence over
the general class default when the Stan code is generated.

## Usage

``` r
bsimms_prior(prior, class = "b", coef = "", resp = "", group = "")
```

## Arguments

- prior:

  Character string: a valid Stan distribution expression, e.g.
  `"normal(0, 1)"`, `"student_t(3, 0, 2.5)"`, `"lkj_corr_cholesky(2)"`,
  or, for `class = "p_global"`, a single positive number (a Dirichlet
  concentration, e.g. `"1"`).

- class:

  One of:

  - `"b"`: fixed-effect slopes (the population-level baseline is not a
    `"b"` coefficient but `"p_global"`, see below).

  - `"p_global"`: Dirichlet concentration for one source's share of the
    global/population-average proportions (MixSIAR's `p_global`).

  - `"sd"`: group-level standard deviations.

  - `"cor"`: group-level correlations (LKJ prior on the Cholesky
    factor).

  - `"sigma"`: residual/observation error; only used for
    `error_structure` `"residual_only"`.

  - `"resid_prop"`: MixSIAR's `resid.prop`, a multiplicative factor
    scaling the propagated source/TDF process variance; only used for
    `error_structure` `"process_residual"`; always restricted to values
    between 0 and 20, matching the range of the recommended default
    prior (see
    [`bsimms_get_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_get_prior.md)),
    so a custom prior can reshape which values in that range are more
    likely but cannot allow values above 20.

  - `"source_mean"`, `"source_sd"`, `"tdf_mean"`, `"tdf_sd"`: only used
    when the corresponding data are supplied raw rather than as
    means/SDs.

  - `"source_cor"`: LKJ prior on the Cholesky factor of each source's
    isotope correlation matrix; only used when source data are raw and
    there are 2+ isotopes.

  - `"resid_cor"`: LKJ prior on the Cholesky factor of the shared
    residual-error correlation matrix; only used for `error_structure`
    `"residual_only"` with 2+ isotopes.

- coef:

  Optional: restrict to one fixed-effect coefficient name, as it appears
  in the model formula's expanded design matrix (e.g. `"SexM"` for a
  factor `Sex`, or `"SexM:RegionB"` for an interaction). Only used with
  `class = "b"`.

- resp:

  Optional: restrict to one isotope (response) name. Used with `class`
  `"sigma"`, `"resid_prop"`, `"source_mean"`, `"source_sd"`,
  `"tdf_mean"` or `"tdf_sd"`.

- group:

  Optional, one of two uses depending on `class`:

  - For `"sd"`/`"cor"`, restrict to one group-level term (the right-hand
    side of a `(... | group)` term, as written in the formula).

  - For `"p_global"`/`"source_mean"`/`"source_sd"`/`"tdf_mean"`/
    `"tdf_sd"`/`"source_cor"`, restrict to one source (as named in
    `source_data`/`tdf_data`), so that source and trophic discrimination
    factor priors can be made source-specific as well as
    isotope-specific. Left unset (`""`), such a prior applies to every
    source (and, except for `"source_cor"`, every isotope).

## Value

A one-row `data.frame` of class `bsimms_prior`.

## Examples

``` r
bsimms_prior("normal(0, 2)", class = "b")
#>         prior class coef resp group
#>  normal(0, 2)     b                
bsimms_prior("normal(1, 1)", class = "b", coef = "SeasonWinter")
#>         prior class         coef resp group
#>  normal(1, 1)     b SeasonWinter           
bsimms_prior("student_t(3, 0, 1)", class = "sd", group = "Region")
#>               prior class coef resp  group
#>  student_t(3, 0, 1)    sd           Region
bsimms_prior(
  "normal(3.4, 0.3)",
  class = "tdf_mean",
  resp = "d15N",
  group = "Beaver"
)
#>             prior    class coef resp  group
#>  normal(3.4, 0.3) tdf_mean      d15N Beaver
bsimms_prior(
  "student_t(3, 0, 0.5)",
  class = "source_sd",
  resp = "d13C",
  group = "Beaver"
)
#>                 prior     class coef resp  group
#>  student_t(3, 0, 0.5) source_sd      d13C Beaver
bsimms_prior("2", class = "p_global", group = "Beaver")
#>  prior    class coef resp  group
#>      2 p_global           Beaver
bsimms_prior("lkj_corr_cholesky(2)", class = "source_cor", group = "Beaver")
#>                 prior      class coef resp  group
#>  lkj_corr_cholesky(2) source_cor           Beaver
bsimms_prior("lkj_corr_cholesky(2)", class = "resid_cor")
#>                 prior     class coef resp group
#>  lkj_corr_cholesky(2) resid_cor                
```
