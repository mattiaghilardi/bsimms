# Specifying source data and error structures

## Overview

[`vignette("getting-started")`](https://mattiaghilardi.github.io/bsimms/articles/getting-started.md)
fits and interprets a first **bsimms** model, using
[`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)’s
default settings with an intercept-only (`~1`) formula. For an
individual mixture sample i and isotope j (superscript m for mixture, s
for source, t for trophic discrimination factor (TDF) throughout), the
expected value is

\mu^m\_{i, j} = \sum\_{k=1}^{K} p\_{i, k} \left(\mu^s\_{k, j} +
\mu^t\_{k, j}\right), \qquad \sum\_{k=1}^{K} p\_{i, k} = 1, \\ \\ p\_{i,
k} \> 0, \tag{1}

where p\_{i, k} is source k’s proportional contribution to sample i, and
\mu^s\_{k, j} and \mu^t\_{k, j} are source k’s mean isotope value and
mean TDF for isotope j. The observed value is normally distributed
around it:

y^m\_{i, j} \sim \mathrm{Normal}\left(\mu^m\_{i, j}, \\ {\sigma^m\_{i,
j}}^2\right). \tag{2}

This vignette works through the modelling options that vary these
equations, beyond those defaults:

- raw vs. summarised source/TDF data
- concentration dependence
- the three error structures
- cross-tracer covariance, automatically estimated where applicable
  rather than a separate option to set

See
[`vignette("covariates")`](https://mattiaghilardi.github.io/bsimms/articles/covariates.md)
for letting source proportions depend on mixture-level covariates,
[`vignette("priors")`](https://mattiaghilardi.github.io/bsimms/articles/priors.md)
for **bsimms**’s data-scaled default priors and how to override them,
and
[`vignette("diagnostics")`](https://mattiaghilardi.github.io/bsimms/articles/diagnostics.md)
for interpreting Stan’s sampler diagnostics.

## Source and TDF data: raw vs. summarised

`source_data` and `tdf_data` can each be supplied in one of two shapes,
controlled by `source_means_sds`/`tdf_means_sds`:

``` r

# raw: one row per replicate measurement, one column per isotope
data.frame(
  Source = rep(c("Beaver", "Deer"), each = 2),
  d13C = c(-25.1, -24.8, -18.2, -17.9),
  d15N = c(5.1, 4.9, 8.2, 7.8)
)
#>   Source  d13C d15N
#> 1 Beaver -25.1  5.1
#> 2 Beaver -24.8  4.9
#> 3   Deer -18.2  8.2
#> 4   Deer -17.9  7.8

# summarised: one row per source, mean/SD columns per isotope
data.frame(
  Source = c("Beaver", "Deer"),
  d13C_mean = c(-25, -18), d13C_sd = c(1, 1),
  d15N_mean = c(5, 8), d15N_sd = c(1, 1)
)
#>   Source d13C_mean d13C_sd d15N_mean d15N_sd
#> 1 Beaver       -25       1         5       1
#> 2   Deer       -18       1         8       1
```

With raw data (`source_means_sds = FALSE`, the default), each replicate
measurement r of source k’s isotope j is hierarchically fitted as

y^s\_{k, j, r} \sim \mathrm{Normal}\left(\mu^s\_{k, j}, \\
{\sigma^s\_{k, j}}^2\right), \tag{3}

with \mu^s\_{k, j} from Equation (1) and \sigma^s\_{k, j} from the
error-structure equations below themselves estimated, so their
uncertainty is fully propagated into the source-proportion posterior.
With 2+ isotopes this is really the diagonal, no-correlation view of a
multivariate model (see [Cross-tracer
covariance](#cross-tracer-covariance) below). With summarised data
(`TRUE`), both are instead fixed at the supplied values. `tdf_data`
works the same way via `tdf_means_sds` (default `TRUE`, since
TDF/fractionation values are usually taken from the literature as fixed
mean/SD pairs, as in MixSIAR) – **bsimms** additionally extends this
hierarchical fitting to raw TDF data, which MixSIAR does not support:

y^t\_{k, j, r} \sim \mathrm{Normal}\left(\mu^t\_{k, j}, \\
{\sigma^t\_{k, j}}^2\right),

with \mu^t\_{k, j} and \sigma^t\_{k, j} likewise estimated when raw,
from Equation (1) and the error-structure equations below respectively.

## Concentration dependence

Concentration dependence ([Phillips & Koch
2002](https://doi.org/10.1007/s004420100786)) weights source proportions
by each source’s elemental concentration for a given tracer, rather than
treating every source as contributing equally regardless of
concentration. It requires an extra `<isotope>_conc` column per isotope
in `source_data` (values in `(0, 1]`, summing to at most 1 per source
across isotopes), and is enabled with `conc_dep = TRUE`. Formally,
p\_{i, k} in Equation (1) is replaced by a concentration-reweighted
p\_{i, k}^\*:

p\_{i, k}^\* = \frac{p\_{i, k} \\ c\_{k, j}}{\sum\_{l=1}^{K} p\_{i, l}
\\ c\_{l, j}},

where c\_{k, j} is source k’s elemental concentration for isotope j (the
`<isotope>_conc` column).

``` r

library(bsimms)

sim_conc <- simulate_bsimms_data(
  ~1,
  n_mixture_obs = 10,
  source_names = c("Beaver", "Deer", "Hare"),
  isotope_names = c("d13C", "d15N"),
  conc_dep = TRUE,
  seed = 1
)

fit_conc <- bsimm(
  sim_conc$formula,
  mixture_data = sim_conc$mixture_data,
  source_data = sim_conc$source_data,
  tdf_data = sim_conc$tdf_data,
  isotope_names = sim_conc$isotope_names,
  source_means_sds = sim_conc$source_means_sds,
  tdf_means_sds = sim_conc$tdf_means_sds,
  conc_dep = sim_conc$conc_dep,
  error_structure = sim_conc$error_structure,
  source_col = sim_conc$source_col,
  chains = 2,
  iter_warmup = 500,
  iter_sampling = 500,
  seed = 1
)
```

``` r

head(sim_conc$source_data)
#>   Source      d13C      d15N d13C_conc d15N_conc
#> 1 Beaver -2.642662 -1.650032 0.3459537  0.423806
#> 2 Beaver -1.381133 -2.701489 0.3459537  0.423806
#> 3 Beaver -1.994469  1.082988 0.3459537  0.423806
#> 4 Beaver -3.277176 -2.558542 0.3459537  0.423806
#> 5 Beaver -2.066832 -1.470974 0.3459537  0.423806
#> 6 Beaver -3.934112 -2.458223 0.3459537  0.423806
fit_conc
#> Bayesian stable isotope mixing model (bsimms)
#>  formula:         ~1
#>  sources (K):     3 (Beaver, Deer, Hare)
#>  isotopes (J):    2 (d13C, d15N)
#>  mixtures (N):    10
#>  error structure: process_residual
#>  source data:     raw
#>  tdf data:        summary
#>  concentration dependence: yes
#>  backend:         cmdstanr
```

## Error structures

`error_structure` controls {\sigma^m\_{i, j}}^2 in Equation (2):

- `"process_residual"` (default): source/TDF variance is propagated into
  the mixture and scaled by an estimated `resid_prop` (\phi_j) factor –
  MixSIAR’s combined `"Residual * Process"` error structure ([Stock &
  Semmens 2016](https://doi.org/10.1002/ecy.1517)): {\sigma^m\_{i, j}}^2
  = \phi_j \sum\_{k=1}^{K} p\_{i, k}^2 \left({\sigma^s\_{k, j}}^2 +
  {\sigma^t\_{k, j}}^2\right).
- `"process_only"`: only source/TDF variance is propagated, with no
  further scaling, as in the original MixSIR model ([Moore & Semmens
  2008](https://doi.org/10.1111/j.1461-0248.2008.01163.x)):
  {\sigma^m\_{i, j}}^2 = \sum\_{k=1}^{K} p\_{i, k}^2
  \left({\sigma^s\_{k, j}}^2 + {\sigma^t\_{k, j}}^2\right).
- `"residual_only"`: source/TDF variance is not propagated into the
  mixture at all; all unexplained variance instead sits directly in
  {\sigma^m_j}^2 – unlike the other two, not indexed by i, since it is
  fit directly rather than derived from p\_{i, k} – an isotope-specific
  `sigma` term, as in the original SIAR model ([Parnell et
  al. 2010](https://doi.org/10.1371/journal.pone.0009672)).

``` r

sim_resid <- simulate_bsimms_data(
  ~1,
  n_mixture_obs = 10,
  source_names = c("Beaver", "Deer", "Hare"),
  isotope_names = c("d13C", "d15N"),
  error_structure = "residual_only",
  seed = 1
)

fit_resid <- bsimm(
  sim_resid$formula,
  mixture_data = sim_resid$mixture_data,
  source_data = sim_resid$source_data,
  tdf_data = sim_resid$tdf_data,
  isotope_names = sim_resid$isotope_names,
  source_means_sds = sim_resid$source_means_sds,
  tdf_means_sds = sim_resid$tdf_means_sds,
  conc_dep = sim_resid$conc_dep,
  error_structure = sim_resid$error_structure,
  source_col = sim_resid$source_col,
  chains = 2,
  iter_warmup = 500,
  iter_sampling = 500,
  seed = 1
)
```

``` r

summary(fit_resid)
#> Bayesian stable isotope mixing model
#>  formula: ~1
#> 
#> Population-average source proportions:
#>  source  mean    sd  q2.5 q97.5  rhat ess_bulk ess_tail
#>  Beaver 0.203 0.094 0.021 0.404 1.023      269      186
#>    Deer 0.399 0.060 0.275 0.508 1.020      330      246
#>    Hare 0.399 0.036 0.324 0.468 1.022      241      187
#> 
#> Error term(s):
#>  isotope  mean    sd  q2.5 q97.5  rhat ess_bulk ess_tail
#>     d13C 0.542 0.139 0.342 0.886 1.003     1446      765
#>     d15N 0.374 0.100 0.239 0.624 1.009      981      522
```

The error term is now `sigma` rather than `resid_prop`.

## Cross-tracer covariance

With 2+ isotopes, cross-tracer covariance ([Hopkins & Ferguson
2012](https://doi.org/10.1371/journal.pone.0028478)) generalises
Equation (2) to vector/matrix form:

\mathbf{y}^m_i \sim
\mathrm{MultivariateNormal}\left(\boldsymbol{\mu}^m_i, \\
\Sigma^m_i\right),

where \Sigma^m_i is the J \times J mixture covariance matrix for sample
i; its diagonal is {\sigma^m\_{i, j}}^2 from the error-structure
equations above. When source data are raw, Equation (3) likewise
generalises:

\mathbf{y}^s\_{k, r} \sim
\mathrm{MultivariateNormal}\left(\boldsymbol{\mu}^s_k, \\
\Sigma^s_k\right),

where \Sigma^s_k is source k’s J \times J covariance matrix; its
diagonal is {\sigma^s\_{k, j}}^2 and its off-diagonal entries are
\sigma^s\_{k, j_1} \sigma^s\_{k, j_2} \\ \rho^s\_{k, j_1, j_2}, with
\rho^s\_{k, j_1, j_2} source k’s estimated correlation between isotopes
j_1 and j_2
([`bsimms_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_prior.md)’s
`"source_cor"` class).

Two specific combinations also estimate \Sigma^m_i’s off-diagonal,
cross-tracer entries automatically:

- Raw source data under a process-based error structure
  (`"process_only"`/`"process_residual"`): source correlation
  (\Sigma^s_k above) propagates into \Sigma^m_i’s off-diagonal the same
  way its diagonal does, \Sigma^m_i\[j_1, j_2\] = \sum\_{k=1}^{K} p\_{i,
  k}^2 \\ \Sigma^s_k\[j_1, j_2\] for j_1 \neq j_2 (also scaled by
  \sqrt{\phi\_{j_1} \phi\_{j_2}} under `"process_residual"`).
- `error_structure = "residual_only"`: a shared residual correlation,
  \Sigma^m\[j_1, j_2\] = \sigma^m\_{j_1} \sigma^m\_{j_2} \\
  \rho^{\text{resid}}\_{j_1, j_2} for j_1 \neq j_2, with \Sigma^m the
  shared J \times J mixture covariance matrix (no i, as for \sigma^m_j
  above) and \rho^{\text{resid}}\_{j_1, j_2} the estimated correlation
  between isotopes j_1 and j_2
  ([`bsimms_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_prior.md)’s
  `"resid_cor"` class).

Neither ever includes a TDF contribution, even when TDF data are raw
(matching MixSIAR, which also never modelled TDF cross-tracer
correlation).

Since source data are raw and there are two isotopes in `fit_resid`
above, its per-source (`Lcorr_source`) and shared residual
(`Lcorr_resid`) correlation Cholesky factors were also estimated; both
are ordinary Stan parameters, retrievable like any other via
[`bsimms_draws()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_draws.md),
but are not part of [`summary()`](https://rdrr.io/r/base/summary.html)’s
output since they are a data/error-structure side effect rather than a
quantity of direct scientific interest here:

``` r

posterior::variables(bsimms_draws(fit_resid, variable = "Lcorr_source"))
#>  [1] "Lcorr_source[1,1,1]" "Lcorr_source[2,1,1]" "Lcorr_source[3,1,1]"
#>  [4] "Lcorr_source[1,2,1]" "Lcorr_source[2,2,1]" "Lcorr_source[3,2,1]"
#>  [7] "Lcorr_source[1,1,2]" "Lcorr_source[2,1,2]" "Lcorr_source[3,1,2]"
#> [10] "Lcorr_source[1,2,2]" "Lcorr_source[2,2,2]" "Lcorr_source[3,2,2]"
posterior::variables(bsimms_draws(fit_resid, variable = "Lcorr_resid"))
#> [1] "Lcorr_resid[1,1]" "Lcorr_resid[2,1]" "Lcorr_resid[1,2]" "Lcorr_resid[2,2]"
```

## References

Phillips DL, Koch PL (2002). “Incorporating concentration dependence in
stable isotope mixing models.” *Oecologia*, 130(1), 114-125.
[doi:10.1007/s004420100786](https://doi.org/10.1007/s004420100786)

Stock BC, Semmens BX (2016). “Unifying error structures in commonly used
biotracer mixing models.” *Ecology*, 97(10), 2952-2960.
[doi:10.1002/ecy.1517](https://doi.org/10.1002/ecy.1517)

Moore JW, Semmens BX (2008). “Incorporating uncertainty and prior
information into stable isotope mixing models.” *Ecology Letters*,
11(5), 470-480.
[doi:10.1111/j.1461-0248.2008.01163.x](https://doi.org/10.1111/j.1461-0248.2008.01163.x)

Parnell AC, Inger R, Bearhop S, Jackson AL (2010). “Source partitioning
using stable isotopes: coping with too much variation.” *PLoS ONE*,
5(3), e9672.
[doi:10.1371/journal.pone.0009672](https://doi.org/10.1371/journal.pone.0009672)

Hopkins JB, Ferguson JM (2012). “Estimating the diets of animals using
stable isotopes and a comprehensive Bayesian mixing model.” *PLoS ONE*,
7(1), e28478.
[doi:10.1371/journal.pone.0028478](https://doi.org/10.1371/journal.pone.0028478)
