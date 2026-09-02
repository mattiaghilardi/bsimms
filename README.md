
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bsimms

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/mattiaghilardi/bsimms/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mattiaghilardi/bsimms/actions/workflows/R-CMD-check.yaml)
[![test_coverage](https://github.com/mattiaghilardi/bsimms/actions/workflows/test_coverage.yaml/badge.svg)](https://github.com/mattiaghilardi/bsimms/actions/workflows/test_coverage.yaml)
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
