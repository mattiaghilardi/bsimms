#' Trace and density plots of model parameters
#'
#' Plots the posterior density and MCMC trace of the model's underlying
#' Stan parameters (population-average source proportions `p_global`,
#' fixed-effect coefficients, group-level standard deviations, and error
#' term(s)), via [bayesplot::mcmc_combo()]. Use [plot_proportions()] for
#' summaries of the source proportions themselves, or
#' [conditional_effects()] to see how they vary with a covariate.
#'
#' @param x A `bsimms_fit` object (as returned by [bsimm()]).
#' @param variable Optional character vector of parameter (base) names to
#'   plot. `NULL` (default) plots `p_global`, the fixed effects (if any),
#'   group-level standard deviations (if any), and the error term(s).
#' @param combo Character vector of two `bayesplot` `mcmc_*` plot types to
#'   combine (default `c("dens", "trace")`); see [bayesplot::mcmc_combo()].
#' @param nvariables Maximum number of parameters shown per plot; models
#'   with more are split across multiple plots (default 5).
#' @param plot Logical; display each plot as a side effect (default
#'   `TRUE`)? If `FALSE`, the plots are only built and returned.
#' @param ask Logical; prompt the user before displaying each new page
#'   after the first (default `TRUE`)? Only relevant if `plot = TRUE` and
#'   there is more than one page.
#' @param newpage Logical; start the first page on a new graphics page
#'   (default `TRUE`)? Every page after the first always starts on a new
#'   page regardless. Only relevant if `plot = TRUE`.
#' @param ... Further arguments passed to [bayesplot::mcmc_combo()].
#' @return A list of plot objects (one per page), invisibly; each is also
#'   displayed as a side effect if `plot = TRUE`.
#' @export
plot.bsimms_fit <- function(x, variable = NULL, combo = c("dens", "trace"), nvariables = 5,
                             plot = TRUE, ask = TRUE, newpage = TRUE, ...) {
  rlang::check_installed(
    "bayesplot",
    reason = "to use `plot.bsimms_fit()` (or use `plot_proportions()` / `summary()` directly)."
  )
  spec <- x$spec
  if (is.null(variable)) {
    variable <- c(
      "p_global",
      if (spec$P > 0) "beta",
      if (length(spec$re_terms) > 0) {
        vapply(spec$re_terms, function(re) paste0("sd_re_", re$label), character(1))
      },
      if (spec$error_structure == "residual_only") "sigma",
      if (spec$error_structure == "process_residual") "resid_prop"
    )
  }

  draws <- bsimms_draws(x, variable = variable)
  all_vars <- posterior::variables(draws)
  if (length(all_vars) == 0) {
    cli::cli_abort("No matching parameters found.", call = NULL)
  }

  if (plot) {
    default_ask <- grDevices::devAskNewPage()
    on.exit(grDevices::devAskNewPage(default_ask))
    grDevices::devAskNewPage(ask = FALSE)
  }

  n_pages <- ceiling(length(all_vars) / nvariables)
  plots <- vector("list", n_pages)
  for (i in seq_len(n_pages)) {
    idx <- ((i - 1) * nvariables + 1):min(i * nvariables, length(all_vars))
    sub_draws <- posterior::subset_draws(draws, variable = all_vars[idx])
    plots[[i]] <- bayesplot::mcmc_combo(sub_draws, combo = combo, ...)
    if (plot) {
      graphics::plot(plots[[i]], newpage = newpage || i > 1)
      if (i == 1) grDevices::devAskNewPage(ask = ask)
    }
  }
  invisible(plots)
}

#' Posterior predictive checks
#'
#' Compares observed mixture isotope values against posterior predictive
#' draws (`y_rep`) using [bayesplot::pp_check()]'s `ppc_*` plotting
#' functions. Since these expect a single response vector, one isotope must
#' be selected.
#'
#' @param object A `bsimms_fit` object (as returned by [bsimm()]).
#' @param resp Character; the isotope to plot (one of `isotope_names`).
#'   Required if the model has more than one isotope; defaults to the only
#'   isotope otherwise.
#' @param type Character; the `bayesplot` PPC plot type, e.g.
#'   `"dens_overlay"` (default), `"hist"`, `"stat"`, `"scatter_avg"`,
#'   `"intervals"` — see [bayesplot::available_ppc()] for the full list
#'   (passed without its `"ppc_"` prefix).
#' @param ndraws Optional integer; number of `y_rep` draws to (randomly)
#'   subsample for the plot. `NULL` (default) uses every draw.
#' @param ... Further arguments passed on to the underlying `ppc_*`
#'   function, e.g. `group` for grouped types.
#' @return A ggplot object, as returned by the underlying `ppc_*` function.
#' @exportS3Method bayesplot::pp_check
pp_check.bsimms_fit <- function(object, resp = NULL, type = "dens_overlay", ndraws = NULL, ...) {
  rlang::check_installed("bayesplot", reason = "to use `pp_check()`.")
  spec <- object$spec
  if (is.null(resp)) {
    if (spec$J > 1) {
      cli::cli_abort(
        "Model has multiple isotopes ({.val {spec$isotope_names}}); specify {.arg resp} to select one.",
        call = NULL
      )
    }
    resp <- spec$isotope_names
  } else {
    resp <- rlang::arg_match0(resp, spec$isotope_names)
  }

  valid_types <- sub("^ppc_", "", as.character(bayesplot::available_ppc("")))
  if (!type %in% valid_types) {
    cli::cli_abort("{.arg type} must be one of {.val {valid_types}}.", call = NULL)
  }
  ppc_fun <- getExportedValue("bayesplot", paste0("ppc_", type))

  j <- match(resp, spec$isotope_names)
  y <- object$standata$y[, j]
  dm <- draws_matrix(object, variable = "y_rep")
  yrep_arr <- extract_array_draws(dm, "y_rep", spec$N, spec$J)
  yrep <- matrix(yrep_arr[, , j], nrow = dim(yrep_arr)[1])

  if (!is.null(ndraws) && ndraws < nrow(yrep)) {
    yrep <- yrep[sample.int(nrow(yrep), ndraws), , drop = FALSE]
  }

  ppc_fun(y = y, yrep = yrep, ...)
}
