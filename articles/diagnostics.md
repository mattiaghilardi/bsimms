# Diagnosing and resolving convergence problems

## Overview

A **bsimms** model is, under the hood, a Stan program fit by Hamiltonian
Monte Carlo (HMC/NUTS), so the same checks that apply to any Stan model
apply here: the sampler’s own warnings (divergent transitions, maximum
treedepth, low E-BFMI) and convergence diagnostics (`rhat`, effective
sample size), on top of the posterior predictive checks already covered
in
[`vignette("getting-started")`](https://mattiaghilardi.github.io/bsimms/articles/getting-started.md).
This vignette works through both, and the causes/fixes most specific to
**bsimms** models.

## Sampler diagnostics

We start by simulating and fitting a simple model:

``` r

library(bsimms)

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
  cores = 4,
  seed = 1
)
```

[`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)
returns the raw backend fit object as `fit$fit`, so the underlying
sampler diagnostics are available exactly as `cmdstanr`/ `rstan` report
them:

``` r

if (fit$backend == "cmdstanr") {
  fit$fit$diagnostic_summary()
} else {
  rstan::check_hmc_diagnostics(fit$fit)
}
#> $num_divergent
#> [1] 0 0 0 0
#> 
#> $num_max_treedepth
#> [1] 0 0 0 0
#> 
#> $ebfmi
#> [1] 0.9778086 0.8496679 0.9392410 0.8887845
```

The three things to look for (see Stan’s own [guide to runtime warnings
and convergence
problems](https://mc-stan.org/learn-stan/diagnostics-warnings.html) for
the full detail):

- **Divergent transitions**: the simulated trajectory the sampler
  integrates diverging from the true, energy-conserving one it
  approximates – usually because the posterior’s curvature varies too
  sharply, in some region, for the step size to resolve. Unlike hitting
  the maximum treedepth, this can bias the resulting draws rather than
  just waste computation: a handful alongside otherwise healthy
  `rhat`/ESS is often tolerable, but more of them, or divergences
  concentrated in a particular region of parameter space, mean that
  region is not being explored correctly and should not be ignored.
- **Maximum treedepth**: the sampler being cut off before it would have
  stopped on its own, to cap how long a single iteration can run. This
  only affects efficiency, not validity – if it is the only warning and
  `rhat`/ESS both look healthy, it is generally safe to ignore.
- **Low E-BFMI**: a validity, not just efficiency, concern – indicating
  that warmup/adaptation did not work well or that the posterior has
  thick tails the sampler is not exploring properly, which can leave
  uncertainty understated in ways `rhat`/effective sample size alone do
  not reveal.

## Convergence: R-hat and effective sample size

`summary(fit)` reports `rhat` and bulk/tail effective sample size (ESS)
alongside the posterior estimates, for the population-average
proportions and any fixed effects, group-level SDs, and error term(s) –
not every underlying Stan parameter. For others (e.g. the raw source/TDF
means/SDs, or the per-observation proportions), pass their name(s) to
`plot(fit, variable = ...)`, or extract them with
[`bsimms_draws()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_draws.md)
and check
[`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html)
directly.

``` r

summary(fit)
#> Bayesian stable isotope mixing model
#>  formula: ~1
#> 
#> Population-average source proportions:
#>  source  mean    sd  q2.5 q97.5  rhat ess_bulk ess_tail
#>  Beaver 0.256 0.120 0.048 0.522 1.003     1485      891
#>    Deer 0.368 0.074 0.201 0.495 1.003     1568     1031
#>    Hare 0.376 0.050 0.270 0.464 1.002     1685     1202
#> 
#> Error term(s):
#>  isotope  mean    sd  q2.5 q97.5  rhat ess_bulk ess_tail
#>     d13C 4.246 2.914 1.125 12.84 1.002     2513     1807
#>     d15N 1.826 1.576 0.477  5.60 1.001     2706     1952
```

As a rule of thumb ([Vehtari et al.
2021](https://doi.org/10.1214/20-BA1221)), treat `rhat > 1.01` or a
bulk/tail ESS below 100 times the number of chains (combined across
chains; 400 for
[`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)’s
default `chains = 4`) as a sign the chains have not mixed well enough to
trust the posterior summaries.

`plot(fit)` visualises the same chains as trace and density plots, which
can reveal *why* (a stuck chain, a multimodal posterior) in a way the
two numbers alone cannot:

``` r

plot(fit)
```

![Density and trace plots of p_global and resid_prop across four MCMC
chains; densities are unimodal and traces overlap with no trend,
indicating good mixing.](diagnostics-plot-fit-1.png)

## Common causes and fixes

- **Sparse group levels with correlated random effects.** A group-level
  term’s slopes/intercepts across ILR dimensions are, by default, fit
  with an estimated correlation. With few levels per group, that
  correlation’s geometry can be poorly identified and cause divergences.
  Raising the target acceptance rate often resolves it – passed via
  [`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)’s
  `...`, as `adapt_delta` for `backend = "cmdstanr"` or
  `control = list(adapt_delta = ...)` for `backend = "rstan"`:

  ``` r

  bsimm(..., adapt_delta = 0.99) # cmdstanr
  bsimm(..., control = list(adapt_delta = 0.99)) # rstan
  ```

  More warmup iterations (`iter_warmup`) can also help the sampler adapt
  its step size to this harder geometry.

- **Frequent maximum-treedepth warnings.** If this is the only warning
  and `rhat`/ESS are otherwise healthy, it is generally safe to leave
  alone. Raising the limit (`max_treedepth` for `cmdstanr`,
  `control = list(max_treedepth = ...)` for `rstan`) trades runtime for
  deeper trajectories, but is not usually the right first response –
  treat frequent hits as a hint to look for the same causes as
  divergences (sparse groups, weak identifiability) rather than raising
  the limit by default.

- **Low effective sample size on the error term.** With little
  replication, decomposing unexplained variance into source/TDF process
  error and observation-level residual error
  (`error_structure = "process_residual"`, the default) can be weakly
  identified. Simplifying to `"residual_only"` (see
  [`vignette("data-error")`](https://mattiaghilardi.github.io/bsimms/articles/data-error.md))
  or gathering more mixture samples both help.

Some restrictions are instead enforced by
[`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)
before fitting even starts, as ordinary errors rather than sampler
diagnostics – for example, `error_structure` being incompatible with a
single mixture data point, or with one mixture sample per level of a
covariate/group. These already explain themselves (each carries an “i”
bullet with the reason), so no separate sampler-side diagnosis is needed
for them.

## References

Vehtari A, Gelman A, Simpson D, Carpenter B, Bürkner PC (2021).
“Rank-normalization, folding, and localization: an improved R-hat for
assessing convergence of MCMC (with discussion).” *Bayesian Analysis*,
16(2), 667-718.
[doi:10.1214/20-BA1221](https://doi.org/10.1214/20-BA1221)
