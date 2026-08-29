#' Bayesian Stable Isotope Mixing Models using Stan
#'
#' `bsimms` fits Bayesian stable isotope mixing models (SIMMs) via
#' [Stan](https://mc-stan.org), estimating the proportional contribution of
#' two or more sources to a mixture from tracer data such as stable
#' isotopes.
#'
#' `bsimms` generates a bespoke Stan program from a user-supplied model
#' specification and fits it via MCMC, using either the `cmdstanr` or
#' `rstan` package as backend. It supports the core modelling options of
#' the JAGS-based package `MixSIAR` (raw or summarised source and trophic
#' discrimination factor data, concentration dependence, flexible error
#' structures, cross-tracer covariance) and extends them by letting source
#' proportions depend on an arbitrary `lme4`-style fixed- and random-effects
#' formula of mixture-level covariates, with weakly informative, data-scaled
#' default priors that can be fully overridden. Source proportions are
#' modelled in isometric log-ratio (ILR) coordinates, so that an
#' unconstrained linear predictor maps onto the source simplex via a
#' numerically stable softmax transform.
#'
#' Main functions:
#' \itemize{
#'   \item [bsimm()] — fit a model (builds + compiles + samples).
#'   \item [make_stancode()] / [make_standata()] — inspect or hand-edit
#'     the generated Stan program / data without fitting.
#'   \item [bsimms_get_prior()] / [bsimms_prior()] — inspect and set priors.
#'   \item [summary.bsimms_fit()], [print.bsimms_fit()] — summarise a
#'     fitted model.
#'   \item [posterior_proportions()], [fitted_proportions()] — posterior
#'     (predictive) source proportions.
#'   \item [posterior_epred.bsimms_fit()], [posterior_predict.bsimms_fit()],
#'     [fitted.bsimms_fit()], [predict.bsimms_fit()] — posterior (predictive)
#'     mixture isotope values.
#'   \item [conditional_effects()], [plot_proportions()],
#'     [plot.bsimms_fit()], [pp_check.bsimms_fit()] — plot a fitted model.
#'   \item [loo.bsimms_fit()], [waic.bsimms_fit()], [loo_compare.bsimms_fit()],
#'     [add_criterion()], [bayes_R2.bsimms_fit()] — model comparison and
#'     evaluation.
#'   \item [draws_long()] — reshape a posterior draws array to long format.
#'   \item [ilr()], [ilr_inv()], [ilr_basis()], [clr()], [clr_inv()] — the
#'     ILR/CLR transforms used throughout, exported for independent use.
#' }
#'
#' @references Carpenter, B., Gelman, A., Hoffman, M.D., Lee, D., Goodrich,
#'   B., Betancourt, M., Brubaker, M., Guo, J., Li, P., & Riddell, A.
#'   (2017). Stan: A probabilistic programming language. *Journal of
#'   Statistical Software*, 76(1), 1-32. \doi{10.18637/jss.v076.i01}
#' @references Stan Development Team. *Stan Modeling Language User's Guide
#'   and Reference Manual*. \url{https://mc-stan.org/docs/}
#' @references Stan Development Team. *RStan: the R interface to Stan*.
#'   \url{https://mc-stan.org/rstan/}
#' @references Gabry, J., Češnovar, R., Johnson, A., & Bronder, S.
#'   *cmdstanr: R Interface to 'CmdStan'*. \url{https://mc-stan.org/cmdstanr/}
#' @references Stock, B.C., Jackson, A.L., Ward, E.J., Parnell, A.C.,
#'   Phillips, D.L., & Semmens, B.X. (2018). Analyzing mixing systems using
#'   a new generation of Bayesian tracer mixing models. *PeerJ*, 6, e5096.
#'   \doi{10.7717/peerj.5096}
#' @keywords internal
#' @importFrom rlang .data
"_PACKAGE"
