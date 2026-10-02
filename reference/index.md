# Package index

## Fit a model

Build, compile, and sample from a bsimms model.

- [`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)
  : Fit a Bayesian stable isotope mixing model

- [`make_stancode()`](https://mattiaghilardi.github.io/bsimms/reference/make_stancode.md)
  :

  Generate the Stan code for a `bsimms` model

- [`make_standata()`](https://mattiaghilardi.github.io/bsimms/reference/make_standata.md)
  :

  Generate the Stan data list for a `bsimms` model

- [`parse_bsimms_formula()`](https://mattiaghilardi.github.io/bsimms/reference/parse_bsimms_formula.md)
  :

  Parse an `lme4`-style formula for the source-proportion model

## Priors

Preview, override, and print bsimms’s data-scaled default priors.

- [`bsimms_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_prior.md)
  :

  Specify a prior for one or more `bsimms` parameters

- [`bsimms_get_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_get_prior.md)
  :

  Default priors for a `bsimms` model

- [`print(`*`<bsimms_prior>`*`)`](https://mattiaghilardi.github.io/bsimms/reference/print.bsimms_prior.md)
  :

  Print a `bsimms_prior` specification

## Inspect a fit

Print and summarise a fitted model.

- [`print(`*`<bsimms_fit>`*`)`](https://mattiaghilardi.github.io/bsimms/reference/print.bsimms_fit.md)
  :

  Print a `bsimms` fit

- [`summary(`*`<bsimms_fit>`*`)`](https://mattiaghilardi.github.io/bsimms/reference/summary.bsimms_fit.md)
  :

  Summarise a `bsimms` fit

## Posterior draws

Extract and reshape raw posterior draws.

- [`bsimms_draws()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_draws.md)
  :

  Extract posterior draws from a `bsimms` fit

- [`draws_long()`](https://mattiaghilardi.github.io/bsimms/reference/draws_long.md)
  : Reshape a draws array into long format

## Predictions

Posterior (predictive) source proportions and mixture isotope values,
for the fitted data or new observations.

- [`posterior_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/posterior_proportions.md)
  : Posterior source proportions
- [`fitted_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/fitted_proportions.md)
  : Posterior source proportions (summarised)
- [`fitted(`*`<bsimms_fit>`*`)`](https://mattiaghilardi.github.io/bsimms/reference/fitted.bsimms_fit.md)
  : Fitted mixture isotope values (summarised)
- [`predict(`*`<bsimms_fit>`*`)`](https://mattiaghilardi.github.io/bsimms/reference/predict.bsimms_fit.md)
  : Posterior predictive draws of mixture isotope values (summarised)
- [`posterior_epred(`*`<bsimms_fit>`*`)`](https://mattiaghilardi.github.io/bsimms/reference/posterior_epred.bsimms_fit.md)
  : Expected mixture isotope values
- [`posterior_predict(`*`<bsimms_fit>`*`)`](https://mattiaghilardi.github.io/bsimms/reference/posterior_predict.bsimms_fit.md)
  : Posterior predictive draws of mixture isotope values

## Model comparison and evaluation

Approximate leave-one-out cross-validation, WAIC, and Bayesian
R-squared.

- [`loo(`*`<bsimms_fit>`*`)`](https://mattiaghilardi.github.io/bsimms/reference/loo.bsimms_fit.md)
  : Approximate leave-one-out cross-validation

- [`loo_compare(`*`<bsimms_fit>`*`)`](https://mattiaghilardi.github.io/bsimms/reference/loo_compare.bsimms_fit.md)
  :

  Compare `bsimms` models by approximate leave-one-out cross-validation

- [`waic(`*`<bsimms_fit>`*`)`](https://mattiaghilardi.github.io/bsimms/reference/waic.bsimms_fit.md)
  : Widely applicable information criterion (WAIC)

- [`add_criterion()`](https://mattiaghilardi.github.io/bsimms/reference/add_criterion.md)
  :

  Cache LOO/WAIC/R-squared criteria on a `bsimms` fit

- [`bayes_R2(`*`<bsimms_fit>`*`)`](https://mattiaghilardi.github.io/bsimms/reference/bayes_R2.bsimms_fit.md)
  : Bayesian R-squared

## Plots

Plot the data before fitting, and diagnostic and summary plots for a
fitted model.

- [`plot(`*`<bsimms_fit>`*`)`](https://mattiaghilardi.github.io/bsimms/reference/plot.bsimms_fit.md)
  : Trace and density plots of model parameters
- [`pp_check(`*`<bsimms_fit>`*`)`](https://mattiaghilardi.github.io/bsimms/reference/pp_check.bsimms_fit.md)
  : Posterior predictive checks
- [`plot_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/plot_proportions.md)
  : Plot posterior source proportions
- [`plot_isospace()`](https://mattiaghilardi.github.io/bsimms/reference/plot_isospace.md)
  : Plot mixture data in isotope space alongside the sources
- [`conditional_effects()`](https://mattiaghilardi.github.io/bsimms/reference/conditional_effects.md)
  [`plot(`*`<bsimms_conditional_effects>`*`)`](https://mattiaghilardi.github.io/bsimms/reference/conditional_effects.md)
  : Conditional effects of covariates on source proportions or predicted
  mixture isotope values

## Simulate data

Simulate data, with known true parameter values, for any bsimms model
configuration.

- [`simulate_bsimms_data()`](https://mattiaghilardi.github.io/bsimms/reference/simulate_bsimms_data.md)
  :

  Simulate data for a `bsimms` model

## ILR/CLR transforms

Transforms between the source simplex and unconstrained coordinates.

- [`ilr()`](https://mattiaghilardi.github.io/bsimms/reference/ilr.md)
  [`ilr_inv()`](https://mattiaghilardi.github.io/bsimms/reference/ilr.md)
  : Forward and inverse ILR transforms
- [`ilr_basis()`](https://mattiaghilardi.github.io/bsimms/reference/ilr_basis.md)
  : Isometric log-ratio (ILR) basis matrix
- [`clr()`](https://mattiaghilardi.github.io/bsimms/reference/clr.md)
  [`clr_inv()`](https://mattiaghilardi.github.io/bsimms/reference/clr.md)
  : Centred log-ratio transform and its inverse

## Package

- [`bsimms`](https://mattiaghilardi.github.io/bsimms/reference/bsimms-package.md)
  [`bsimms-package`](https://mattiaghilardi.github.io/bsimms/reference/bsimms-package.md)
  : Bayesian Stable Isotope Mixing Models using Stan
