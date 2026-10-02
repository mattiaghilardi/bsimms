# Bayesian Stable Isotope Mixing Models using Stan

`bsimms` fits Bayesian stable isotope mixing models (SIMMs) via
[Stan](https://mc-stan.org), estimating the proportional contribution of
two or more sources to a mixture from tracer data such as stable
isotopes.

## Details

`bsimms` generates a bespoke Stan program from a user-supplied model
specification and fits it via MCMC, using either the `cmdstanr` or
`rstan` package as backend. It supports the core modelling options of
the JAGS-based package `MixSIAR` (raw or summarised source and trophic
discrimination factor data, concentration dependence, flexible error
structures, cross-tracer covariance) and extends them by letting source
proportions depend on an arbitrary `lme4`-style fixed- and
random-effects formula of mixture-level covariates, with weakly
informative, data-scaled default priors that can be fully overridden.
Source proportions are modelled in isometric log-ratio (ILR)
coordinates, so that an unconstrained linear predictor maps onto the
source simplex via a numerically stable softmax transform.

Main functions:

- [`bsimm()`](https://mattiaghilardi.github.io/bsimms/reference/bsimm.md)
  — fit a model (builds + compiles + samples).

- [`make_stancode()`](https://mattiaghilardi.github.io/bsimms/reference/make_stancode.md)
  /
  [`make_standata()`](https://mattiaghilardi.github.io/bsimms/reference/make_standata.md)
  — inspect or hand-edit the generated Stan program / data without
  fitting.

- [`bsimms_get_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_get_prior.md)
  /
  [`bsimms_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_prior.md)
  — inspect and set priors.

- [`plot_isospace()`](https://mattiaghilardi.github.io/bsimms/reference/plot_isospace.md)
  — plot the mixtures against the sources in isotope space, before
  fitting.

- [`summary.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/summary.bsimms_fit.md),
  [`print.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/print.bsimms_fit.md)
  — summarise a fitted model.

- [`posterior_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/posterior_proportions.md),
  [`fitted_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/fitted_proportions.md)
  — posterior (predictive) source proportions.

- [`posterior_epred.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/posterior_epred.bsimms_fit.md),
  [`posterior_predict.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/posterior_predict.bsimms_fit.md),
  [`fitted.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/fitted.bsimms_fit.md),
  [`predict.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/predict.bsimms_fit.md)
  — posterior (predictive) mixture isotope values.

- [`conditional_effects()`](https://mattiaghilardi.github.io/bsimms/reference/conditional_effects.md),
  [`plot_proportions()`](https://mattiaghilardi.github.io/bsimms/reference/plot_proportions.md),
  [`plot.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/plot.bsimms_fit.md),
  [`pp_check.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/pp_check.bsimms_fit.md)
  — plot a fitted model.

- [`loo.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/loo.bsimms_fit.md),
  [`waic.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/waic.bsimms_fit.md),
  [`loo_compare.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/loo_compare.bsimms_fit.md),
  [`add_criterion()`](https://mattiaghilardi.github.io/bsimms/reference/add_criterion.md),
  [`bayes_R2.bsimms_fit()`](https://mattiaghilardi.github.io/bsimms/reference/bayes_R2.bsimms_fit.md)
  — model comparison and evaluation.

- [`draws_long()`](https://mattiaghilardi.github.io/bsimms/reference/draws_long.md)
  — reshape a posterior draws array to long format.

- [`simulate_bsimms_data()`](https://mattiaghilardi.github.io/bsimms/reference/simulate_bsimms_data.md)
  — simulate data for any model configuration, with known true parameter
  values.

- [`ilr()`](https://mattiaghilardi.github.io/bsimms/reference/ilr.md),
  [`ilr_inv()`](https://mattiaghilardi.github.io/bsimms/reference/ilr.md),
  [`ilr_basis()`](https://mattiaghilardi.github.io/bsimms/reference/ilr_basis.md),
  [`clr()`](https://mattiaghilardi.github.io/bsimms/reference/clr.md),
  [`clr_inv()`](https://mattiaghilardi.github.io/bsimms/reference/clr.md)
  — the ILR/CLR transforms used throughout, exported for independent
  use.

## References

Carpenter, B., Gelman, A., Hoffman, M.D., Lee, D., Goodrich, B.,
Betancourt, M., Brubaker, M., Guo, J., Li, P., & Riddell, A. (2017).
Stan: A probabilistic programming language. *Journal of Statistical
Software*, 76(1), 1-32.
[doi:10.18637/jss.v076.i01](https://doi.org/10.18637/jss.v076.i01)

Stan Development Team. *Stan Modeling Language User's Guide and
Reference Manual*. <https://mc-stan.org/docs/>

Stan Development Team. *RStan: the R interface to Stan*.
<https://mc-stan.org/rstan/>

Gabry, J., Češnovar, R., Johnson, A., & Bronder, S. *cmdstanr: R
Interface to 'CmdStan'*. <https://mc-stan.org/cmdstanr/>

Stock, B.C., Jackson, A.L., Ward, E.J., Parnell, A.C., Phillips, D.L., &
Semmens, B.X. (2018). Analyzing mixing systems using a new generation of
Bayesian tracer mixing models. *PeerJ*, 6, e5096.
[doi:10.7717/peerj.5096](https://doi.org/10.7717/peerj.5096)

## See also

Useful links:

- <https://github.com/mattiaghilardi/bsimms>

- Report bugs at <https://github.com/mattiaghilardi/bsimms/issues>

## Author

**Maintainer**: Mattia Ghilardi <mattia.ghilardi91@gmail.com>
([ORCID](https://orcid.org/0000-0001-9592-7252))

Authors:

- Mattia Ghilardi <mattia.ghilardi91@gmail.com>
  ([ORCID](https://orcid.org/0000-0001-9592-7252))
