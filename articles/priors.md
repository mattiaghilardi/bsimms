# Setting custom priors

## Overview

[`vignette("data-error")`](https://mattiaghilardi.github.io/bsimms/articles/data-error.md)
and
[`vignette("covariates")`](https://mattiaghilardi.github.io/bsimms/articles/covariates.md)
work through **bsimms**’s modelling options. This vignette covers its
priors. Every parameter a
[`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)
model has – source proportions, fixed- and random-effect coefficients,
error terms, and, when source/TDF data are raw, their means/SDs – gets a
default prior automatically. Each can be previewed with
[`bsimms_get_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_get_prior.md)
before fitting, and overridden, in full or in part, with
[`bsimms_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_prior.md).

Some of these defaults are scaled to the data instead of fixed: the
residual error term under `error_structure = "residual_only"` (`sigma`),
and, when source/TDF data are raw, each source’s own mean and SD
(`source_mean`/`source_sd`, `tdf_mean`/`tdf_sd`). These use the sample
median and MAD (median absolute deviation) of the relevant data, rather
than a single fixed default shared across every model and dataset.
Median and MAD are robust to the skew and outliers common in isotope
data.

## Previewing default priors

[`bsimms_get_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_get_prior.md)
previews the default priors **bsimms** would use for a given model
configuration, without fitting it:

``` r

library(bsimms)

sim <- simulate_bsimms_data(
  ~1,
  n_mixture_obs = 10,
  source_names = c("Beaver", "Deer", "Hare"),
  isotope_names = c("d13C", "d15N"),
  seed = 1
)

bsimms_get_prior(
  ~1,
  mixture_data = sim$mixture_data,
  source_data = sim$source_data,
  tdf_data = sim$tdf_data,
  isotope_names = sim$isotope_names,
  error_structure = sim$error_structure
)
#>                      prior       class coef resp  group
#>               normal(0, 1)           b                 
#>                          1    p_global           Beaver
#>                          1    p_global             Deer
#>                          1    p_global             Hare
#>             uniform(0, 20)  resid_prop                 
#>             uniform(0, 20)  resid_prop      d13C       
#>             uniform(0, 20)  resid_prop      d15N       
#>  normal(-2.33985, 9.42643) source_mean      d13C Beaver
#>  student_t(3, 0, 0.942643)   source_sd      d13C Beaver
#>  normal(-7.32937, 10.5985) source_mean      d13C   Deer
#>   student_t(3, 0, 1.05985)   source_sd      d13C   Deer
#>   normal(12.8653, 9.38194) source_mean      d13C   Hare
#>  student_t(3, 0, 0.938194)   source_sd      d13C   Hare
#>   normal(-2.34629, 13.109) source_mean      d15N Beaver
#>    student_t(3, 0, 1.3109)   source_sd      d15N Beaver
#>  normal(-12.6946, 6.53279) source_mean      d15N   Deer
#>  student_t(3, 0, 0.653279)   source_sd      d15N   Deer
#>   normal(7.55435, 16.8137) source_mean      d15N   Hare
#>   student_t(3, 0, 1.68137)   source_sd      d15N   Hare
#>       lkj_corr_cholesky(1)  source_cor                 
#>       lkj_corr_cholesky(1)  source_cor           Beaver
#>       lkj_corr_cholesky(1)  source_cor             Deer
#>       lkj_corr_cholesky(1)  source_cor             Hare
```

## Overriding priors

Any row can be overridden with
[`bsimms_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_prior.md),
passed to
[`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)’s
`prior` argument. A more specific row (naming a `coef`/`resp`/`group`)
takes precedence over a general `class` default, and a general override
cascades onto every more specific default row it would otherwise be
shadowed by, rather than being silently unused.
[`bsimms_get_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_get_prior.md)
above pre-fills one `p_global` row per source; a general override with
no `group` applies to all three at once:

``` r

fit_general <- bsimm(
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
  prior = bsimms_prior("2", class = "p_global"),
  chains = 2,
  iter_warmup = 500,
  iter_sampling = 500,
  seed = 1
)
```

``` r

fit_general$prior[fit_general$prior$class == "p_global", ]
#>  prior    class coef resp  group
#>      2 p_global           Beaver
#>      2 p_global             Deer
#>      2 p_global             Hare
#>      2 p_global
```

A `group`-specific override instead only replaces that one source’s row,
taking precedence over the general one. For example, to additionally
give Beaver’s share a stronger-than-default concentration:

``` r

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
  prior = bsimms_prior("3", class = "p_global", group = "Beaver"),
  chains = 2,
  iter_warmup = 500,
  iter_sampling = 500,
  seed = 1
)
```

``` r

fit$prior[fit$prior$class == "p_global" & fit$prior$group == "Beaver", ]
#>  prior    class coef resp  group
#>      3 p_global           Beaver
```

Several overrides can be combined with
[`c()`](https://rdrr.io/r/base/c.html), e.g. to also raise Deer’s
concentration at the same time:

``` r

c(
  bsimms_prior("3", class = "p_global", group = "Beaver"),
  bsimms_prior("2", class = "p_global", group = "Deer")
)
#>  prior    class coef resp  group
#>      3 p_global           Beaver
#>      2 p_global             Deer
```

Other parameter classes work the same way.
[`bsimms_get_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_get_prior.md)’s
output above shows `resid_prop`’s default is `uniform(0, 20)` for both
isotopes; `resid_prop` is restricted to non-negative values, so a custom
prior needs a matching, positive-only distribution – e.g. a
`gamma(2, 2)`, centred near its expected value of 1, for `d13C`:

``` r

fit_resid_prop <- bsimm(
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
  prior = bsimms_prior("gamma(2, 2)", class = "resid_prop", resp = "d13C"),
  chains = 2,
  iter_warmup = 500,
  iter_sampling = 500,
  seed = 1
)
```

``` r

fit_resid_prop$prior[fit_resid_prop$prior$class == "resid_prop", ]
#>           prior      class coef resp group
#>  uniform(0, 20) resid_prop                
#>     gamma(2, 2) resid_prop      d13C      
#>  uniform(0, 20) resid_prop      d15N
```

### Informative priors from independent source-contribution data

When independent information on each source’s contribution already
exists – from a previous study, a different data type, or expert
knowledge – it can be turned into an informative Dirichlet prior for
every source’s `p_global` at once. Rescale the raw counts or percentages
so they sum to K, the number of sources ([Stock et al.
2018](https://doi.org/10.7717/peerj.5096), eq. 4): this keeps the
prior’s total weight the same as the default, “generalist”
\mathrm{Dirichlet}(1, \ldots, 1) prior, while preserving the independent
data’s relative proportions among sources.

A common example is diet composition estimated from stomach or faecal
content analysis, but the same approach applies to any
source-contribution estimate from outside the isotope data itself –
e.g. known upstream contributions in a sediment- or pollutant-tracing
application:

``` r

# hypothetical stomach-content counts for each source
diet_counts <- c(Beaver = 30, Deer = 18, Hare = 12)

# rescale so the concentrations sum to K (Stock et al. 2018, eq. 4)
alpha <- length(diet_counts) * diet_counts / sum(diet_counts)
alpha
#> Beaver   Deer   Hare 
#>    1.5    0.9    0.6
```

``` r

informative_prior <- c(
  bsimms_prior("1.5", class = "p_global", group = "Beaver"),
  bsimms_prior("0.9", class = "p_global", group = "Deer"),
  bsimms_prior("0.6", class = "p_global", group = "Hare")
)

fit_informative <- bsimm(
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
  prior = informative_prior,
  chains = 2,
  iter_warmup = 500,
  iter_sampling = 500,
  seed = 1
)
```

``` r

fit_informative$prior[fit_informative$prior$class == "p_global", ]
#>  prior    class coef resp  group
#>    1.5 p_global           Beaver
#>    0.9 p_global             Deer
#>    0.6 p_global             Hare
```

## References

Stock BC, Jackson AL, Ward EJ, Parnell AC, Phillips DL, Semmens BX
(2018). “Analyzing mixing systems using a new generation of Bayesian
tracer mixing models.” *PeerJ*, 6, e5096.
[doi:10.7717/peerj.5096](https://doi.org/10.7717/peerj.5096)
