
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bsimms

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/mattiaghilardi/bsimms/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mattiaghilardi/bsimms/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/mattiaghilardi/bsimms/graph/badge.svg)](https://app.codecov.io/gh/mattiaghilardi/bsimms)
<!-- [![CRAN status](https://www.r-pkg.org/badges/version/bsimms)](https://CRAN.R-project.org/package=bsimms) -->
<!-- badges: end -->

## Overview

bsimms is an R package for fitting Bayesian stable isotope mixing models
(SIMMs) via Stan, estimating the proportional contribution of two or
more sources to a mixture from tracer data such as stable isotopes. It
brings together mixing-model advances developed across the literature
and implemented in the JAGS-based package MixSIAR — hierarchical source
fitting, fractionation error, concentration dependence, cross-tracer
covariance, flexible error structures — and extends them by letting
source proportions depend on an arbitrary lme4-style fixed- and
random-effects formula of mixture-level covariates.

## Installation

You can install the development version of bsimms from GitHub like so:

``` r
# install.packages("pak")
pak::pak("mattiaghilardi/bsimms")
```

Fitting models also requires a working Stan backend, either
[cmdstanr](https://mc-stan.org/cmdstanr/) (recommended) or
[rstan](https://mc-stan.org/rstan/), neither of which is installed
automatically:

``` r
# cmdstanr (not on CRAN)
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::check_cmdstan_toolchain()
cmdstanr::install_cmdstan()

# or rstan (on CRAN)
Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1) # only necessary for Linux without the nodejs library / headers
install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
```

See the [CmdStanR installation
guide](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) or the
[RStan Getting Started
guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) for
OS-specific instructions and troubleshooting.

## How to use bsimms

As an example, we fit a mixing model to simulated data for a 3-source,
2-isotope model with 10 mixture samples, no covariates (`~1`), raw
replicate source data, and summarised (mean/SD) TDF data.

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

Printing the fit gives a quick overview of the model.

``` r
fit
#> Bayesian stable isotope mixing model (bsimms)
#>  formula:         ~1
#>  sources (K):     3 (Beaver, Deer, Hare)
#>  isotopes (J):    2 (d13C, d15N)
#>  mixtures (N):    10
#>  error structure: process_residual
#>  source data:     raw
#>  tdf data:        summary
#>  backend:         cmdstanr
```

`summary()` reports the population-average source proportions and error
term(s), each with a posterior mean, SD, 95% credible interval, and MCMC
convergence diagnostics (`rhat`, bulk/tail effective sample size). With
the default `error_structure = "process_residual"`, the error term shown
per isotope is `resid_prop`, MixSIAR’s multiplicative factor scaling the
source/TDF variance propagated into the mixture – not a raw residual SD.

``` r
summary(fit)
#> Bayesian stable isotope mixing model
#>  formula: ~1
#> 
#> Population-average source proportions:
#>  source  mean    sd  q2.5 q97.5 rhat ess_bulk ess_tail
#>  Beaver 0.256 0.122 0.037 0.513    1     2248     1901
#>    Deer 0.368 0.075 0.207 0.498    1     2200     1783
#>    Hare 0.376 0.051 0.271 0.472    1     2604     2208
#> 
#> Error term(s):
#>  isotope  mean    sd  q2.5  q97.5 rhat ess_bulk ess_tail
#>     d13C 4.189 2.807 1.112 12.150    1     4127     2847
#>     d15N 1.787 1.410 0.469  5.656    1     3420     2414
```

`rhat` close to 1 and large effective sample sizes indicate the chains
mixed well, which we can also check visually with trace and density
plots of the underlying Stan parameters.

``` r
plot(fit)
```

<img src="man/figures/README-plot-fit-1.png" alt="Density and trace plots of p_global and resid_prop across four MCMC chains; densities are unimodal and traces overlap with no trend, indicating good mixing." width="70%" style="display: block; margin: auto;" />

Population-average contributions lean towards Deer and Hare over Beaver,
though with fairly wide credible intervals given the small simulated
sample size.

Since this model has no covariates, every mixture sample shares the same
population-average proportions, so we can visualise their posterior
distribution directly.

``` r
p_arr <- posterior_proportions(fit)
plot_proportions(p_arr[, 1, , drop = FALSE], type = "density")
```

<img src="man/figures/README-plot-proportions-1.png" alt="Density plot of posterior source proportions for Beaver, Deer, and Hare; Beaver's distribution is shifted towards lower proportions and more widely spread than the overlapping Deer and Hare distributions." width="70%" style="display: block; margin: auto;" />
