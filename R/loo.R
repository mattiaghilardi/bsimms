#' Approximate leave-one-out cross-validation
#'
#' Computes PSIS-LOO (Vehtari, Gelman, and Gabry 2017) from the model's
#' `log_lik` generated quantity, which is the joint log density of each
#' mixture sample's full isotope profile (isotopes are summed/modelled
#' jointly per mixture sample, never held out individually, so the mixture
#' sample is the leave-one-out unit -- see [make_stancode()]).
#'
#' @param x A `bsimms_fit` object (as returned by [bsimm()]).
#' @param cores Number of cores used for [loo::relative_eff()] and
#'   [loo::loo()] (default `getOption("mc.cores", 1)`).
#' @param ... Further arguments passed on to [loo::loo.array()].
#' @return A `loo` object, as returned by [loo::loo()].
#' @references Vehtari, A., Gelman, A., & Gabry, J. (2017). Practical
#'   Bayesian model evaluation using leave-one-out cross-validation and WAIC.
#'   *Statistics and Computing*, 27(5), 1413-1432.
#'   \doi{10.1007/s11222-016-9696-4}
#' @exportS3Method loo::loo
#' @examples
#' \donttest{
#' sim <- simulate_bsimms_data(
#'   ~1,
#'   n_mixture_obs = 10,
#'   source_names = c("Beaver", "Deer", "Hare"),
#'   isotope_names = c("d13C", "d15N"),
#'   seed = 1
#' )
#' fit <- bsimm(
#'   sim$formula, mixture_data = sim$mixture_data,
#'   source_data = sim$source_data, tdf_data = sim$tdf_data,
#'   isotope_names = sim$isotope_names,
#'   source_means_sds = sim$source_means_sds, tdf_means_sds = sim$tdf_means_sds,
#'   conc_dep = sim$conc_dep, error_structure = sim$error_structure,
#'   source_col = sim$source_col,
#'   chains = 2, iter_warmup = 500, iter_sampling = 500
#' )
#' loo::loo(fit)
#' }
loo.bsimms_fit <- function(x, cores = getOption("mc.cores", 1), ...) {
  rlang::check_installed("loo", reason = "to use `loo()`.")
  ll <- unclass(bsimms_draws(x, variable = "log_lik")) # iteration x chain x N
  r_eff <- loo::relative_eff(exp(ll), cores = cores)
  loo::loo(ll, r_eff = r_eff, cores = cores, ...)
}

#' Widely applicable information criterion (WAIC)
#'
#' Computes WAIC (Watanabe 2010) from the model's `log_lik` generated
#' quantity (see [loo.bsimms_fit()] for what constitutes a leave-one-out
#' unit here). [loo.bsimms_fit()] is recommended over WAIC since PSIS-LOO
#' additionally provides Pareto k reliability diagnostics; WAIC is provided
#' mainly for comparison with WAIC values reported by other mixing-model
#' software.
#'
#' @param x A `bsimms_fit` object (as returned by [bsimm()]).
#' @param ... Further arguments passed on to [loo::waic.array()].
#' @return A `waic` object, as returned by [loo::waic()].
#' @references Watanabe, S. (2010). Asymptotic equivalence of Bayes cross
#'   validation and widely applicable information criterion in singular
#'   learning theory. *Journal of Machine Learning Research*, 11, 3571-3594.
#' @exportS3Method loo::waic
#' @examples
#' \donttest{
#' sim <- simulate_bsimms_data(
#'   ~1,
#'   n_mixture_obs = 10,
#'   source_names = c("Beaver", "Deer", "Hare"),
#'   isotope_names = c("d13C", "d15N"),
#'   seed = 1
#' )
#' fit <- bsimm(
#'   sim$formula, mixture_data = sim$mixture_data,
#'   source_data = sim$source_data, tdf_data = sim$tdf_data,
#'   isotope_names = sim$isotope_names,
#'   source_means_sds = sim$source_means_sds, tdf_means_sds = sim$tdf_means_sds,
#'   conc_dep = sim$conc_dep, error_structure = sim$error_structure,
#'   source_col = sim$source_col,
#'   chains = 2, iter_warmup = 500, iter_sampling = 500
#' )
#' loo::waic(fit)
#' }
waic.bsimms_fit <- function(x, ...) {
  rlang::check_installed("loo", reason = "to use `waic()`.")
  ll <- unclass(bsimms_draws(x, variable = "log_lik")) # iteration x chain x N
  loo::waic(ll, ...)
}

#' Cache LOO/WAIC/R-squared criteria on a `bsimms` fit
#'
#' Computes one or more model-evaluation criteria and stores them in
#' `x$criteria` (a named list), returning the augmented fit.
#' [loo_compare.bsimms_fit()] reuses a cached `"loo"`/`"waic"` criterion
#' automatically instead of recomputing it; `"bayes_R2"` is cached purely
#' for reuse/comparison.
#'
#' @param x A `bsimms_fit` object (as returned by [bsimm()]).
#' @param criterion Character vector of criteria to compute and cache:
#'   `"loo"` (default), `"waic"`, `"bayes_R2"`, or any combination.
#' @param ... Further arguments passed on to [loo.bsimms_fit()] (if
#'   `"loo"` is requested), [waic.bsimms_fit()] (`"waic"`), or
#'   [bayes_R2.bsimms_fit()] (`"bayes_R2"`, always cached with
#'   `summary = FALSE`, i.e. the raw posterior draws).
#' @return `x`, with `x$criteria` updated to include the newly computed
#'   criterion/criteria.
#' @export
#' @examples
#' \donttest{
#' sim <- simulate_bsimms_data(
#'   ~1,
#'   n_mixture_obs = 10,
#'   source_names = c("Beaver", "Deer", "Hare"),
#'   isotope_names = c("d13C", "d15N"),
#'   seed = 1
#' )
#' fit <- bsimm(
#'   sim$formula, mixture_data = sim$mixture_data,
#'   source_data = sim$source_data, tdf_data = sim$tdf_data,
#'   isotope_names = sim$isotope_names,
#'   source_means_sds = sim$source_means_sds, tdf_means_sds = sim$tdf_means_sds,
#'   conc_dep = sim$conc_dep, error_structure = sim$error_structure,
#'   source_col = sim$source_col,
#'   chains = 2, iter_warmup = 500, iter_sampling = 500
#' )
#' fit <- add_criterion(fit, "loo")
#' fit$criteria$loo
#' }
add_criterion <- function(x, ...) {
  UseMethod("add_criterion")
}

#' @rdname add_criterion
#' @export
add_criterion.bsimms_fit <- function(x, criterion = "loo", ...) {
  criterion <- rlang::arg_match(criterion, c("loo", "waic", "bayes_R2"), multiple = TRUE)
  for (crit in criterion) {
    x$criteria[[crit]] <- switch(crit,
      loo = loo.bsimms_fit(x, ...),
      waic = waic.bsimms_fit(x, ...),
      bayes_R2 = bayes_R2.bsimms_fit(x, summary = FALSE, ...)
    )
  }
  x
}

#' Compare `bsimms` models by approximate leave-one-out cross-validation
#'
#' Computes [loo::loo()] or [loo::waic()] for each model (a `loo`/`waic`
#' object passed directly, e.g. one already computed via `loo(fit)`, is used
#' as-is; for a `bsimms_fit`, a criterion already cached via
#' [add_criterion()] is reused instead of recomputed) and ranks them with
#' [loo::loo_compare()].
#'
#' @param x A `bsimms_fit` object. Must be a `bsimms_fit` (not a
#'   precomputed `loo`/`waic` object) for this method to be dispatched to;
#'   a precomputed `loo`/`waic` object can only appear in `...`.
#' @param ... Additional `bsimms_fit` objects and/or `loo`/`waic` objects to
#'   compare against `x`.
#' @param criterion `"loo"` (default) or `"waic"`; which criterion to use
#'   for any `bsimms_fit` objects passed in (cached via [add_criterion()]
#'   if present, computed on the fly otherwise). Ignored for arguments that
#'   are already `loo`/`waic` objects.
#' @param model_names Optional character vector naming the compared models,
#'   in the order `x`, `...`. Defaults to the deparsed argument expressions.
#' @return The result of [loo::loo_compare()]: model names are given in its
#'   `model` column, ranked by expected log pointwise predictive density
#'   (best model first).
#' @exportS3Method loo::loo_compare
#' @examples
#' \donttest{
#' sim <- simulate_bsimms_data(
#'   ~Sex,
#'   n_mixture_obs = 10,
#'   source_names = c("Beaver", "Deer", "Hare"),
#'   isotope_names = c("d13C", "d15N"),
#'   n_levels = list(Sex = 2),
#'   seed = 1
#' )
#' fit1 <- bsimm(
#'   ~1, mixture_data = sim$mixture_data,
#'   source_data = sim$source_data, tdf_data = sim$tdf_data,
#'   isotope_names = sim$isotope_names,
#'   source_means_sds = sim$source_means_sds, tdf_means_sds = sim$tdf_means_sds,
#'   conc_dep = sim$conc_dep, error_structure = sim$error_structure,
#'   source_col = sim$source_col,
#'   chains = 2, iter_warmup = 500, iter_sampling = 500
#' )
#' fit2 <- bsimm(
#'   ~Sex, mixture_data = sim$mixture_data,
#'   source_data = sim$source_data, tdf_data = sim$tdf_data,
#'   isotope_names = sim$isotope_names,
#'   source_means_sds = sim$source_means_sds, tdf_means_sds = sim$tdf_means_sds,
#'   conc_dep = sim$conc_dep, error_structure = sim$error_structure,
#'   source_col = sim$source_col,
#'   chains = 2, iter_warmup = 500, iter_sampling = 500
#' )
#' loo::loo_compare(fit1, fit2)
#' }
loo_compare.bsimms_fit <- function(x, ..., criterion = c("loo", "waic"), model_names = NULL) {
  criterion <- rlang::arg_match(criterion)
  rlang::check_installed("loo", reason = "to use `loo_compare()`.")
  models <- c(list(x), list(...))
  if (is.null(model_names)) {
    call_args <- as.list(match.call())[-1]
    call_args <- call_args[names(call_args) %in% c("", "x")]
    model_names <- vapply(call_args, rlang::expr_deparse, character(1))
  }

  matches_criterion <- if (criterion == "waic") loo::is.waic else loo::is.psis_loo
  precomputed <- !vapply(models, inherits, logical(1), "bsimms_fit")
  mismatched <- precomputed & !vapply(models, matches_criterion, logical(1))
  if (any(mismatched)) {
    cli::cli_abort(
      "Cannot compare models evaluated with different criteria: {.val {model_names[mismatched]}} must be {.cls {criterion}} object{?s} to match {.code criterion = '{criterion}'}.",
      call = NULL
    )
  }

  compute <- if (criterion == "loo") loo.bsimms_fit else waic.bsimms_fit
  crits <- Map(function(m, nm) {
    if (!inherits(m, "bsimms_fit")) return(m)
    cached <- m$criteria[[criterion]]
    if (!is.null(cached)) return(cached)
    cli::cli_inform(
      "Computing {.val {criterion}} for model {.val {nm}} (not cached; use {.fn add_criterion} to cache it for reuse)."
    )
    compute(m)
  }, models, model_names)
  names(crits) <- model_names
  loo::loo_compare(crits)
}
