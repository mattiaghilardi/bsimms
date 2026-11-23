# Modelling covariate effects

## Overview

[`vignette("getting-started")`](https://mattiaghilardi.github.io/bsimms/articles/getting-started.md)
and
[`vignette("data-error")`](https://mattiaghilardi.github.io/bsimms/articles/data-error.md)
fit **bsimms** models with an intercept-only (`~1`) `formula`. With no
covariates, p\_{i, k} in Equation (1) of
[`vignette("data-error")`](https://mattiaghilardi.github.io/bsimms/articles/data-error.md)
is the same for every sample i: a single, population-average value,
`p_global`, drawn from a Dirichlet distribution. This vignette covers
`formula` itself: an `lme4`-style specification of mixture-level
covariates, spanning fixed effects (continuous, factor, interaction) and
crossed/nested random effects, that instead let p\_{i, k} vary across
samples.

Covariate effects are added, in isometric log-ratio (ILR) coordinates,
as deviations from this `p_global` baseline. `p_global` is first
forward-transformed into ILR coordinates, then combined with the
covariate terms and mapped back onto the source simplex:

\begin{aligned} \eta\_{\text{global}} &=
\mathrm{ilr}(p\_{\text{global}}), \\ \eta_i &= \eta\_{\text{global}} +
\mathbf{x}\_i^\top \boldsymbol{\beta} + \mathbf{z}\_i^\top
\mathbf{b}\_{g(i)}, \\ p_i &= \mathrm{ilr}^{-1}(\eta_i). \end{aligned}

where \eta_i and \eta\_{\text{global}} are each vectors of D = K - 1 ILR
coordinates, with \eta\_{\text{global}} playing the role of the
intercept in this per-dimension linear model since neither \mathbf{x}\_i
(sample i’s row of fixed-effect covariates) nor \boldsymbol{\beta}
include an intercept term of their own; \mathbf{z}\_i and
\mathbf{b}\_{g(i)} are the analogous random-effect design row and
group-level deviation for the group g(i) sample i belongs to (one such
term per grouping factor, when `formula` has more than one); and
\mathrm{ilr}^{-1} maps \eta_i back onto the source simplex, giving p_i,
a vector of K source proportions – the p\_{i, k} used in Equation (1) of
[`vignette("data-error")`](https://mattiaghilardi.github.io/bsimms/articles/data-error.md).

By default (`center = TRUE`, `scale = TRUE` in
[`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)/
[`make_stancode()`](https://mattiaghilardi.github.io/bsimms/reference/make_stancode.md)/[`make_standata()`](https://mattiaghilardi.github.io/bsimms/reference/make_standata.md)),
every continuous fixed-effect column of \mathbf{x}\_i is centred and
scaled before fitting – factor/dummy-coded columns are left untouched.
This means `p_global` represents the population-average proportions when
every continuous covariate is at its mean value, not at a raw value of
zero, which can be a nonsensical point to extrapolate to (e.g. a
consumer’s body length of zero).
[`fitted_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/fitted_proportions.md)/[`predict()`](https://rdrr.io/r/stats/predict.html)/[`conditional_effects()`](https://mattiaghilardi.github.io/bsimms/reference/conditional_effects.md)
all report values on the original covariate scale regardless, so this
only affects how `p_global` and \boldsymbol{\beta} themselves are
interpreted, not predictions at any particular covariate value.

## Fixed and random effects

Here we simulate a factor fixed effect (`Sex`), a continuous fixed
effect (`Age`), their interaction, and a random intercept (`Region`, 3
levels):

``` r

library(bsimms)

sim_cov <- simulate_bsimms_data(
  ~ Sex * Age + (1 | Region),
  n_mixture_obs = 60,
  source_names = c("Beaver", "Deer", "Hare"),
  isotope_names = c("d13C", "d15N"),
  n_levels = list(Sex = 2),
  n_groups = list(Region = 3),
  seed = 1
)

fit_full <- bsimm(
  sim_cov$formula,
  mixture_data = sim_cov$mixture_data,
  source_data = sim_cov$source_data,
  tdf_data = sim_cov$tdf_data,
  isotope_names = sim_cov$isotope_names,
  source_means_sds = sim_cov$source_means_sds,
  tdf_means_sds = sim_cov$tdf_means_sds,
  conc_dep = sim_cov$conc_dep,
  error_structure = sim_cov$error_structure,
  source_col = sim_cov$source_col,
  cores = 4,
  seed = 1
)
```

[`print()`](https://rdrr.io/r/base/print.html) confirms that `Age`, the
only continuous fixed-effect covariate here, was centred and scaled by
default (see the Overview above):

``` r

fit_full
#> Bayesian stable isotope mixing model (bsimms)
#>  formula:         ~Sex * Age + (1 | Region)
#>  sources (K):     3 (Beaver, Deer, Hare)
#>  isotopes (J):    2 (d13C, d15N)
#>  mixtures (N):    60
#>  error structure: process_residual
#>  source data:     raw
#>  tdf data:        summary
#>  covariates:      centred and scaled
#>  group-level terms:
#>    ((Intercept) | Region), 3 levels
#>  backend:         cmdstanr
```

[`summary()`](https://rdrr.io/r/base/summary.html) now also reports
fixed-effect coefficients (ILR-scale deviations from `p_global`) and
group-level standard deviations:

``` r

summary(fit_full)
#> Bayesian stable isotope mixing model
#>  formula: ~Sex * Age + (1 | Region)
#> 
#> Population-average source proportions:
#>  source  mean    sd  q2.5 q97.5  rhat ess_bulk ess_tail
#>  Beaver 0.338 0.094 0.161 0.542 1.004     1501     1573
#>    Deer 0.384 0.095 0.196 0.579 1.003     1592     1662
#>    Hare 0.278 0.078 0.143 0.452 1.005     1179     1550
#> 
#> Fixed effects (ILR scale, deviations from the global proportions):
#>  coefficient ilr_dim   mean    sd   q2.5  q97.5  rhat ess_bulk ess_tail
#>         SexB       1  0.976 0.245  0.536  1.501 1.001     2221     2088
#>          Age       1  0.003 0.096 -0.187  0.193 1.000     1914     2269
#>     SexB:Age       1 -0.639 0.175 -0.989 -0.300 1.001     2223     2338
#>         SexB       2 -0.732 0.112 -0.958 -0.516 1.001     2027     2420
#>          Age       2  0.158 0.056  0.044  0.267 1.001     2011     1891
#>     SexB:Age       2  0.057 0.105 -0.146  0.269 1.000     1785     2115
#> 
#> Group-level standard deviations (ILR scale):
#>   group        term ilr_dim  mean    sd  q2.5 q97.5  rhat ess_bulk ess_tail
#>  Region (Intercept)       1 0.567 0.380 0.138 1.563 1.002     1503     1616
#>  Region (Intercept)       2 0.527 0.314 0.160 1.317 1.001     1369     1388
#> 
#> Error term(s):
#>  isotope  mean    sd  q2.5 q97.5  rhat ess_bulk ess_tail
#>     d13C 2.929 0.925 1.484 5.032 1.001     1983     1911
#>     d15N 1.897 0.553 1.017 3.167 1.002     2383     2100
```

## Nested and crossed random effects

Nested random effects (e.g. sampling sites within regions) use `lme4`’s
`/` shorthand: `(1 | Region/Site)` expands to
`(1 | Region) + (1 | Site:Region)`, giving one random intercept per
`Region` and a further one per `Region`-`Site` combination. Crossed
random effects, where the grouping factors have no parent/child
relationship (e.g. `(1 | Sex:Region)`), use `:` the same way, without
the implied nesting. `chains`/`iter_warmup`/`iter_sampling` are reduced
below purely to keep this demonstration fast.

``` r

sim_nested <- simulate_bsimms_data(
  ~ (1 | Region / Site),
  n_mixture_obs = 60,
  source_names = c("Beaver", "Deer", "Hare"),
  isotope_names = c("d13C", "d15N"),
  n_groups = list(Region = 3, Site = 6),
  seed = 1
)

fit_nested <- bsimm(
  sim_nested$formula,
  mixture_data = sim_nested$mixture_data,
  source_data = sim_nested$source_data,
  tdf_data = sim_nested$tdf_data,
  isotope_names = sim_nested$isotope_names,
  source_means_sds = sim_nested$source_means_sds,
  tdf_means_sds = sim_nested$tdf_means_sds,
  conc_dep = sim_nested$conc_dep,
  error_structure = sim_nested$error_structure,
  source_col = sim_nested$source_col,
  chains = 2,
  iter_warmup = 500,
  iter_sampling = 500,
  seed = 1
)
```

[`summary()`](https://rdrr.io/r/base/summary.html) now reports a
separate group-level SD for `Region` and for `Site:Region` – each
region’s own baseline deviation, and each site’s further deviation
within its region. `n_groups$Site`’s total (`6`) is split across
`Region`’s levels the same way `n_mixture_obs` is split across a
factor’s levels elsewhere, here 2 sites per region:

``` r

summary(fit_nested)
#> Bayesian stable isotope mixing model
#>  formula: ~(1 | Region/Site)
#> 
#> Population-average source proportions:
#>  source  mean    sd  q2.5 q97.5  rhat ess_bulk ess_tail
#>  Beaver 0.422 0.127 0.185 0.680 1.003      367      410
#>    Deer 0.395 0.094 0.197 0.583 1.005      438      625
#>    Hare 0.183 0.083 0.039 0.363 1.003      485      356
#> 
#> Group-level standard deviations (ILR scale):
#>        group        term ilr_dim  mean    sd  q2.5 q97.5  rhat ess_bulk
#>  Site:Region (Intercept)       1 0.229 0.206 0.006 0.794 1.003      347
#>  Site:Region (Intercept)       2 0.299 0.335 0.009 1.153 1.001      464
#>       Region (Intercept)       1 0.568 0.433 0.027 1.614 1.003      401
#>       Region (Intercept)       2 0.525 0.505 0.014 1.910 1.003      560
#>  ess_tail
#>       399
#>       644
#>       307
#>       507
#> 
#> Error term(s):
#>  isotope  mean    sd  q2.5 q97.5  rhat ess_bulk ess_tail
#>     d13C 1.726 0.582 0.814 3.122 1.001      763      601
#>     d15N 2.783 0.798 1.506 4.590 1.003      781      794
```

## Model comparison

To check whether `fit_full`’s covariates (`Sex`, `Age`, and `Region`),
taken together, are worth including, fit the same data with the null
(intercept-only) model and compare with `loo_compare()`:

``` r

fit_null_cov <- bsimm(
  ~1,
  mixture_data = sim_cov$mixture_data,
  source_data = sim_cov$source_data,
  tdf_data = sim_cov$tdf_data,
  isotope_names = sim_cov$isotope_names,
  source_means_sds = sim_cov$source_means_sds,
  tdf_means_sds = sim_cov$tdf_means_sds,
  conc_dep = sim_cov$conc_dep,
  error_structure = sim_cov$error_structure,
  source_col = sim_cov$source_col,
  cores = 4,
  seed = 1
)
```

``` r

loo::loo_compare(fit_full, fit_null_cov)
#> Computing "loo" for model "fit_full" (not cached; use `add_criterion()` to
#> cache it for reuse).
#> Computing "loo" for model "fit_null_cov" (not cached; use `add_criterion()` to
#> cache it for reuse).
#>         model elpd_diff se_diff p_worse diag_diff diag_elpd
#>      fit_full       0.0     0.0      NA                    
#>  fit_null_cov     -67.2     5.8    1.00   N < 100
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
```

The best-fitting model is listed first (`elpd_diff = 0`); a positive,
several-times-larger-than-its-SE `elpd_diff` for the runner-up would
indicate a real difference in predictive fit.

## Bayesian R²

With covariates that create real across-sample variation,
[`rstantools::bayes_R2()`](https://mc-stan.org/rstantools/reference/bayes_R2.html)
is now informative (unlike the intercept-only case in
[`vignette("getting-started")`](https://mattiaghilardi.github.io/bsimms/articles/getting-started.md),
where it is always ~0):

``` r

rstantools::bayes_R2(fit_full)
#>   isotope      mean         sd       q2.5     q97.5
#> 1    d13C 0.1805656 0.05989236 0.06898159 0.2991976
#> 2    d15N 0.8929404 0.01127905 0.86540764 0.9083137
```

R² is reported per isotope, and relates directly to that isotope’s error
term in `summary(fit_full)` above: `resid_prop` scales how much
source/TDF process variance is propagated into the mixture as residual
noise, so a smaller `resid_prop` – more of that isotope’s variance
explained by the covariates rather than left as noise – corresponds to a
higher R², and vice versa, as seen here for `d15N` (smaller
`resid_prop`, higher R²) versus `d13C` (larger `resid_prop`, lower R²).

## Conditional effects

[`conditional_effects()`](https://mattiaghilardi.github.io/bsimms/reference/conditional_effects.md)
shows how source proportions vary across a fixed-effect covariate,
holding every other covariate at a reference value and group-level terms
at their population-average level. With no `effects` argument, it
computes one plot per term of the fixed-effect formula – main effects
and two-way interactions alike:

``` r

ce <- conditional_effects(fit_full)
plot(ce)
```

![Forest plot of posterior source proportions for Beaver, Deer, and Hare
across the two levels of Sex, each with 50% and 95% credible
intervals.](covariates-conditional-effects-1.png)![Ribbon plot of
posterior source proportions for Beaver, Deer, and Hare across the range
of Age, each with 50% and 95% credible
intervals.](covariates-conditional-effects-2.png)![Forest plot of
posterior source proportions for Beaver, Deer, and Hare across the two
levels of Sex, faceted by three representative values of Age, each with
50% and 95% credible intervals.](covariates-conditional-effects-3.png)

A factor covariate (`Sex`) is shown as a forest-plot-style linerange at
its levels; a numeric covariate (`Age`) as a ribbon across its observed
range. The interaction (`Sex:Age`) facets `Sex`’s linerange plot by
three representative values of `Age` (its mean and +/-1 SD by default,
overridable via `int_conditions`).

An `effects` entry names which covariate is swept (the first) and which
becomes the moderator (the second), so naming the same interaction the
other way round – `"Age:Sex"` instead of `"Sex:Age"` – instead sweeps
`Age` as a ribbon, faceted by `Sex`’s two levels:

``` r

ce_age <- conditional_effects(fit_full, effects = "Age:Sex")
plot(ce_age)
```

![Ribbon plot of posterior source proportions for Beaver, Deer, and Hare
across the range of Age, faceted by the two levels of Sex, each with 50%
and 95% credible intervals.](covariates-conditional-effects-swap-1.png)

The random intercept on `Region` is, by default, marginalised out of
these plots (`re_formula = NA`); pass a reduced formula
(e.g. `~ (1 | Region)`) to condition on it instead.

## Visualising the random effect

[`conditional_effects()`](https://mattiaghilardi.github.io/bsimms/reference/conditional_effects.md)
only ever sweeps a fixed-effect covariate, so it cannot show `Region`’s
own random effect. To visualise it instead, predict directly for each of
`Region`’s three fitted levels with
[`posterior_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/posterior_proportions.md)/[`fitted_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/fitted_proportions.md)
and `newdata`, holding `Sex`/`Age` at a reference value and conditioning
on `Region` via `re_formula`:

``` r

region_data <- data.frame(
  Sex = levels(sim_cov$mixture_data$Sex)[1],
  Age = mean(sim_cov$mixture_data$Age),
  Region = levels(sim_cov$mixture_data$Region)
)
region_data
#>   Sex         Age Region
#> 1   A -0.04012237      A
#> 2   A -0.04012237      B
#> 3   A -0.04012237      C
```

``` r

p_region <- posterior_proportions(
  fit_full,
  newdata = region_data,
  re_formula = ~ (1 | Region)
)
```

[`plot_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/plot_proportions.md)
labels the x axis generically (“Observation”), since the
`[n_draws, n_obs, K]` array returned by
[`posterior_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/posterior_proportions.md)
carries no information about what varies between rows. It returns a
plain `ggplot` object, though, so relabelling it with `Region`’s own
values needs only a couple of `ggplot2` calls appended:

``` r

plot_proportions(p_region, type = "interval") +
  ggplot2::labs(x = "Region") +
  ggplot2::scale_x_discrete(labels = region_data$Region)
```

![Forest plot of posterior source proportions for Beaver, Deer, and Hare
across the three levels of Region, each with 50% and 95% credible
intervals.](covariates-plot-region-1.png)

## Predicting for unobserved data

The same `newdata` mechanism, introduced in
[`vignette("getting-started")`](https://mattiaghilardi.github.io/bsimms/articles/getting-started.md),
also predicts for entirely new mixture samples – including, here, a
specific individual (`Sex = "B"`, `Age = 2`) from a genuinely new,
previously unobserved `Region`. Set `allow_new_levels = TRUE` to sample
a fresh group-level deviation for the new region instead of erroring;
with as few levels as here, `sample_new_levels = "gaussian"` (drawing
from the estimated group-level SD) is more appropriate than the default
`"uncertainty"` (resampling one of the existing levels’ draws):

``` r

new_region_data <- data.frame(Sex = "B", Age = 2, Region = "D")
p_new_region <- posterior_proportions(
  fit_full,
  newdata = new_region_data,
  re_formula = ~ (1 | Region),
  allow_new_levels = TRUE,
  sample_new_levels = "gaussian"
)
```

``` r

plot_proportions(p_new_region, type = "density")
```

![Density plot of posterior source proportions for Beaver, Deer, and
Hare for a new, previously unobserved
Region.](covariates-plot-new-region-1.png)

The same `newdata`/`re_formula`/`allow_new_levels`/`sample_new_levels`
arguments work identically for
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html)/[`predict()`](https://rdrr.io/r/stats/predict.html),
returning this new individual’s predicted mixture isotope values instead
of its source proportions. Once that individual is actually measured,
this prediction is a benchmark to flag a surprising sample: compare its
observed isotope values against it, the same idea as
[`vignette("getting-started")`](https://mattiaghilardi.github.io/bsimms/articles/getting-started.md)’s
posterior predictive checks, but for one specific new sample rather than
the whole fitted dataset:

``` r

predict(
  fit_full,
  newdata = new_region_data,
  re_formula = ~ (1 | Region),
  allow_new_levels = TRUE,
  sample_new_levels = "gaussian"
)
#>   row isotope     mean       sd       q2.5     q97.5
#> 1   1    d13C 4.144067 2.050610 -0.1768738  8.061565
#> 2   1    d15N 5.734920 3.178086 -0.8252824 12.202831
```
