#' Fit a Bayesian stable isotope mixing model
#'
#' The main entry point. Builds a bespoke Stan program for the requested
#' model (see [make_stancode()]), compiles it, and draws posterior samples
#' via MCMC (`cmdstanr` preferred, `rstan` as a fallback).
#'
#' @inheritParams make_stancode
#' @param backend `"cmdstanr"` or `"rstan"`. Default: `cmdstanr` if
#'   installed, otherwise `rstan`.
#' @param chains Positive integer, number of MCMC chains (default 4).
#' @param iter_warmup Positive integer, warmup iterations per chain (default
#'   1000).
#' @param iter_sampling Positive integer, post-warmup sampling iterations
#'   per chain (default 1000).
#' @param seed Optional integer seed.
#' @param cores Positive integer, number of cores for parallel chains
#'   (default `getOption("mc.cores", 1)`).
#' @param refresh How often to print sampler progress (in iterations);
#'   default lets the backend choose.
#' @param ... Further arguments passed on to `cmdstanr`'s
#'   `$sample()` or `rstan::sampling()`.
#' @return An object of class `bsimms_fit`: a list with elements `fit` (the
#'   raw backend fit object), `backend`, `spec` (internal model
#'   specification), `stancode`, `standata`, `prior`, and `call`.
#' @inherit make_stancode references
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
#' }
bsimm <- function(formula, mixture_data, source_data, tdf_data, isotope_names,
                    source_means_sds = FALSE, tdf_means_sds = TRUE,
                    conc_dep = FALSE,
                    error_structure = c("process_residual", "process_only", "residual_only"),
                    prior = NULL,
                    source_col = "Source",
                    backend = c("auto", "cmdstanr", "rstan"),
                    chains = 4, iter_warmup = 1000, iter_sampling = 1000,
                    seed = NULL, cores = getOption("mc.cores", 1), refresh = NULL,
                    ...) {
  error_structure <- rlang::arg_match(error_structure)
  backend <- rlang::arg_match(backend)
  mc <- match.call()

  spec <- build_bsimms_spec(
    formula = formula, mixture_data = mixture_data,
    source_data = source_data, tdf_data = tdf_data, isotope_names = isotope_names,
    source_means_sds = source_means_sds, tdf_means_sds = tdf_means_sds,
    conc_dep = conc_dep, error_structure = error_structure, source_col = source_col
  )
  prior_df <- merge_bsimms_prior(default_bsimms_prior(spec), prior)
  code <- bsimms_stancode_from_spec(spec, prior_df)
  class(code) <- c("bsimms_stancode", "character")
  sdata <- standata_from_spec(spec, prior_df)
  class(sdata) <- c("bsimms_standata", "list")

  backend <- resolve_backend(backend)

  fit <- if (backend == "cmdstanr") {
    fit_cmdstanr(code, sdata, chains, iter_warmup, iter_sampling, seed, cores, refresh, ...)
  } else {
    fit_rstan(code, sdata, chains, iter_warmup, iter_sampling, seed, cores, refresh, ...)
  }

  structure(
    list(
      fit = fit, backend = backend, spec = spec,
      stancode = code, standata = sdata, prior = prior_df, call = mc
    ),
    class = "bsimms_fit"
  )
}

#' Resolve `bsimm()`'s `backend = "auto"` into a concrete, installed
#' backend (`cmdstanr` preferred, `rstan` as a fallback), or validate that
#' an explicitly requested backend is actually installed.
#'
#' @param backend `"auto"`, `"cmdstanr"`, or `"rstan"`.
#' @return `"cmdstanr"` or `"rstan"`.
#' @noRd
resolve_backend <- function(backend) {
  if (backend == "auto") {
    if (requireNamespace("cmdstanr", quietly = TRUE)) return("cmdstanr")
    if (requireNamespace("rstan", quietly = TRUE)) return("rstan")
    cli::cli_abort(
      c(
        "Fitting requires either {.pkg cmdstanr} or {.pkg rstan} to be installed.",
        "*" = "{.code install.packages(\"rstan\")}",
        "*" = "or",
        "*" = "{.code install.packages(\"cmdstanr\", repos = c(\"https://mc-stan.org/r-packages/\", getOption(\"repos\")))}",
        "*" = "{.code cmdstanr::install_cmdstan()}",
        "i" = "You can still call {.fn make_stancode}/{.fn make_standata} and fit the model yourself."
      ),
      call = NULL
    )
  }
  rlang::check_installed(backend, reason = paste0("to use `backend = \"", backend, "\"`."))
  backend
}

#' Compile and sample a model via `cmdstanr`.
#'
#' @param code A `bsimms_stancode` string (as returned by
#'   `bsimms_stancode_from_spec()`).
#' @param sdata A `bsimms_standata` list (as returned by
#'   `standata_from_spec()`).
#' @param chains,iter_warmup,iter_sampling,seed,cores,refresh Same meaning
#'   as the corresponding arguments of `bsimm()`.
#' @param ... Further arguments passed on to `cmdstanr`'s `$sample()`.
#' @return The `CmdStanMCMC` fit object returned by `$sample()`.
#' @noRd
fit_cmdstanr <- function(code, sdata, chains, iter_warmup, iter_sampling, seed, cores, refresh, ...) {
  stan_file <- cmdstanr::write_stan_file(unclass(code))
  mod <- cmdstanr::cmdstan_model(stan_file)
  args <- list(
    data = unclass(sdata), chains = chains,
    parallel_chains = max(1, cores),
    iter_warmup = iter_warmup, iter_sampling = iter_sampling
  )
  if (!is.null(seed)) args$seed <- seed
  if (!is.null(refresh)) args$refresh <- refresh
  do.call(mod$sample, c(args, list(...)))
}

#' Compile and sample a model via `rstan`.
#'
#' @param code A `bsimms_stancode` string (as returned by
#'   `bsimms_stancode_from_spec()`).
#' @param sdata A `bsimms_standata` list (as returned by
#'   `standata_from_spec()`).
#' @param chains,iter_warmup,iter_sampling,seed,cores,refresh Same meaning
#'   as the corresponding arguments of `bsimm()` (`iter_warmup`/
#'   `iter_sampling` combined into `rstan::stan()`'s `iter`/`warmup`).
#' @param ... Further arguments passed on to `rstan::stan()`.
#' @return The `stanfit` object returned by `rstan::stan()`.
#' @noRd
fit_rstan <- function(code, sdata, chains, iter_warmup, iter_sampling, seed, cores, refresh, ...) {
  args <- list(
    model_code = unclass(code), data = unclass(sdata), chains = chains,
    iter = iter_warmup + iter_sampling, warmup = iter_warmup,
    cores = max(1, cores)
  )
  if (!is.null(seed)) args$seed <- seed
  if (!is.null(refresh)) args$refresh <- refresh
  do.call(rstan::stan, c(args, list(...)))
}
