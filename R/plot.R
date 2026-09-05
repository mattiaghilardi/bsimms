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
#' plot(fit)
#' }
plot.bsimms_fit <- function(
  x,
  variable = NULL,
  combo = c("dens", "trace"),
  nvariables = 5,
  plot = TRUE,
  ask = TRUE,
  newpage = TRUE,
  ...
) {
  rlang::check_installed(
    "bayesplot",
    reason = paste0(
      "to use `plot.bsimms_fit()` (or use `plot_proportions()` / ",
      "`summary()` directly)."
    )
  )
  spec <- x$spec
  if (is.null(variable)) {
    variable <- c(
      "p_global",
      if (spec$P > 0) "beta",
      if (length(spec$re_terms) > 0) {
        vapply(
          spec$re_terms,
          function(re) paste0("sd_re_", re$label),
          character(1)
        )
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
#'   subsample for the plot. `NULL` (default) uses every draw for PPC types
#'   that aggregate/summarise across draws (e.g. `"stat"`, `"intervals"`,
#'   `"scatter_avg"`), where more draws only improve precision, or 10 draws
#'   for types that overlay one `y_rep` dataset per draw (e.g.
#'   `"dens_overlay"`, `"hist"`), where using every draw would overplot the
#'   figure.
#' @param ... Further arguments passed on to the underlying `ppc_*`
#'   function, e.g. `group` for grouped types.
#' @return A ggplot object, as returned by the underlying `ppc_*` function.
#' @exportS3Method bayesplot::pp_check
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
#' bayesplot::pp_check(fit, resp = "d13C")
#' }
pp_check.bsimms_fit <- function(
  object,
  resp = NULL,
  type = "dens_overlay",
  ndraws = NULL,
  ...
) {
  rlang::check_installed("bayesplot", reason = "to use `pp_check()`.")
  spec <- object$spec
  if (is.null(resp)) {
    if (spec$J > 1) {
      cli::cli_abort(
        paste0(
          "Model has multiple isotopes ({.val {spec$isotope_names}}); ",
          "specify {.arg resp} to select one."
        ),
        call = NULL
      )
    }
    resp <- spec$isotope_names
  } else {
    resp <- rlang::arg_match0(resp, spec$isotope_names)
  }

  valid_types <- sub("^ppc_", "", as.character(bayesplot::available_ppc("")))
  if (!type %in% valid_types) {
    cli::cli_abort(
      "{.arg type} must be one of {.val {valid_types}}.",
      call = NULL
    )
  }
  ppc_fun <- getExportedValue("bayesplot", paste0("ppc_", type))

  if (is.null(ndraws)) {
    if (type %in% aggregate_ppc_types) {
      cli::cli_inform(
        "Using all posterior draws for ppc type {.val {type}} by default."
      )
    } else {
      ndraws <- 10
      cli::cli_inform(
        "Using {ndraws} posterior draws for ppc type {.val {type}} by default."
      )
    }
  }

  j <- match(resp, spec$isotope_names)
  y <- object$standata$y[, j]
  dm <- draws_matrix(object, variable = "y_rep")
  dm <- subset_ndraws(dm, ndraws)
  yrep_arr <- extract_array_draws(dm, "y_rep", spec$N, spec$J)
  yrep <- matrix(yrep_arr[,, j], nrow = dim(yrep_arr)[1])

  ppc_fun(y = y, yrep = yrep, ...)
}

#' `bayesplot` PPC types (`ppc_*`, prefix stripped) that aggregate/
#' summarise across draws rather than overlaying one `y_rep` dataset per
#' draw, so [pp_check.bsimms_fit()] can default `ndraws` to "use every
#' draw" only for these.
#' @noRd
aggregate_ppc_types <- c(
  "error_scatter_avg",
  "error_scatter_avg_vs_x",
  "intervals",
  "intervals_grouped",
  "loo_intervals",
  "loo_pit",
  "loo_pit_overlay",
  "loo_pit_qq",
  "loo_ribbon",
  "loo_pit_ecdf",
  "pit_ecdf",
  "pit_ecdf_grouped",
  "ribbon",
  "ribbon_grouped",
  "rootogram",
  "scatter_avg",
  "scatter_avg_grouped",
  "stat",
  "stat_2d",
  "stat_freqpoly_grouped",
  "stat_grouped",
  "violin_grouped"
)

#' Plot posterior source proportions
#'
#' Plots posterior source proportions, as returned by
#' [posterior_proportions()] or [fitted_proportions()] (`summary =
#' FALSE`): a density or histogram of the posterior distribution
#' (`type = "density"`/`"histogram"`, one observation only), or nested
#' credible intervals across one or more observations (`type =
#' "interval"`, a forest/caterpillar plot coloured by source). Requires
#' `ggplot2`.
#'
#' @param p_arr A numeric `[n_draws, n_obs, K]` array of posterior
#'   proportion draws, with source names attached as the 3rd dimension's
#'   `dimnames`, as returned by [posterior_proportions()].
#' @param type Character; `"density"` (default), `"histogram"`, or
#'   `"interval"`. `"density"`/`"histogram"` require `p_arr` to have
#'   exactly one observation (row); use `"interval"` for more than one.
#' @param probs One or more credible-interval masses to display when
#'   `type = "interval"`, e.g. the default `c(0.5, 0.95)` draws both a 50%
#'   and a 95% interval. Ignored for `"density"`/`"histogram"`.
#' @param robust Logical; if `FALSE` (default) the point estimate (for
#'   `type = "interval"`) is the `mean`, if `TRUE` the `median`. Ignored
#'   for `"density"`/`"histogram"`.
#' @param point_size Size of the point marking the central estimate, for
#'   `type = "interval"` (default 2). Ignored otherwise.
#' @param ... Further arguments passed to the underlying `ggplot2` geom:
#'   [ggplot2::geom_density()], [ggplot2::geom_histogram()], or
#'   [ggplot2::geom_linerange()] for `type` `"density"`, `"histogram"`, or
#'   `"interval"` respectively.
#' @return A `ggplot` object.
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
#' p_arr <- posterior_proportions(fit)
#' plot_proportions(p_arr, type = "interval")
#' plot_proportions(p_arr[, 1, , drop = FALSE], type = "density")
#' }
plot_proportions <- function(
  p_arr,
  type = c("density", "histogram", "interval"),
  probs = c(0.5, 0.95),
  robust = FALSE,
  point_size = 2,
  ...
) {
  type <- rlang::arg_match(type)
  rlang::check_installed("ggplot2", reason = "to use `plot_proportions()`.")

  if (length(dim(p_arr)) != 3 || is.null(dimnames(p_arr)[[3]])) {
    cli::cli_abort(
      c(
        paste0(
          "{.arg p_arr} must be a `[n_draws, n_obs, K]` array with ",
          "source names attached as the 3rd dimension's dimnames."
        ),
        "i" = paste0(
          "Use {.fn posterior_proportions} (or {.fn fitted_proportions} ",
          "with {.code summary = FALSE}) to build it."
        )
      ),
      call = NULL
    )
  }

  n_obs <- dim(p_arr)[2]
  if (type %in% c("density", "histogram") && n_obs != 1) {
    cli::cli_abort(
      c(
        paste0(
          "{.arg p_arr} must have exactly one observation (row) for ",
          "{.val {type}} plots, not {n_obs}."
        ),
        "i" = paste0(
          "Use {.code type = \"interval\"} to plot proportions for ",
          "multiple observations."
        )
      ),
      call = NULL
    )
  }

  if (type %in% c("density", "histogram")) {
    plot_proportions_dist(p_arr, type, ...)
  } else {
    plot_proportions_interval(p_arr, probs, robust, point_size, ...)
  }
}

#' Build the `type = "density"`/`"histogram"` plot for `plot_proportions()`:
#' posterior proportions on the x axis, coloured/filled by source.
#'
#' @param p_arr Numeric `[n_draws, 1, K]` array of proportion draws
#'   (`plot_proportions()` requires exactly one observation for these
#'   types), with source names as the 3rd dimension's `dimnames`.
#' @param type `"density"` or `"histogram"`.
#' @param ... Further arguments passed to [ggplot2::geom_density()] or
#'   [ggplot2::geom_histogram()].
#' @return A `ggplot` object.
#' @noRd
plot_proportions_dist <- function(p_arr, type, ...) {
  df <- draws_long(p_arr, var_col = "source", value_col = "proportion")
  ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data$proportion,
      color = .data$source,
      fill = .data$source
    )
  ) +
    (if (type == "density") {
      ggplot2::geom_density(alpha = 0.3, ...)
    } else {
      ggplot2::geom_histogram(alpha = 0.5, position = "identity", ...)
    }) +
    ggplot2::labs(
      x = "Posterior proportion",
      y = if (type == "density") "Density" else "Count",
      color = "Source",
      fill = "Source"
    ) +
    ggplot2::theme_minimal()
}

#' Validate a vector of credible-interval masses, sorted widest-first.
#'
#' @param probs Numeric vector of one or more credible-interval masses;
#'   each must be strictly between 0 and 1 (e.g. `0.95` for a 95%
#'   interval).
#' @return `probs`, sorted decreasing (widest interval first).
#' @noRd
validate_ci_probs <- function(probs) {
  if (length(probs) < 1 || any(probs <= 0 | probs >= 1)) {
    cli::cli_abort(
      paste0(
        "{.arg probs} must be one or more credible-interval masses ",
        "strictly between 0 and 1."
      ),
      call = NULL
    )
  }
  sort(probs, decreasing = TRUE)
}

#' Long-format multi-interval summary of an `[n_draws, n_obs, n_var]`
#' draws array (e.g. source-proportion or isotope-value draws): one row
#' per (observation, category, interval width), with columns `row`,
#' `<cat_col>`, `estimate`, `lower`, `upper`, `width` (an ordered factor,
#' widest first).
#'
#' @param arr Numeric `[n_draws, n_obs, n_var]` array of draws, with
#'   category names (e.g. source or isotope names) as the 3rd dimension's
#'   `dimnames`.
#' @param probs One or more credible-interval masses (validated and
#'   sorted via `validate_ci_probs()`); each produces one `lower`/`upper`
#'   pair per (observation, category).
#' @param robust Logical; if `FALSE` (default) `estimate` is the `mean`,
#'   if `TRUE` the `median`.
#' @param cat_col Name to give the category-identity column (default
#'   `"source"`).
#' @return A long-format data frame; see above for columns.
#' @noRd
summarise_multi_interval <- function(
  arr,
  probs,
  robust = FALSE,
  cat_col = "source"
) {
  probs <- validate_ci_probs(probs)
  cat_names <- dimnames(arr)[[3]]
  est_measure <- if (robust) "median" else "mean"
  n_draws <- dim(arr)[1]
  n_obs <- dim(arr)[2]

  out <- do.call(
    rbind,
    lapply(seq_len(n_obs), function(i) {
      d <- posterior::as_draws_matrix(matrix(
        arr[, i, ],
        nrow = n_draws,
        dimnames = list(NULL, cat_names)
      ))
      est <- as.data.frame(posterior::summarise_draws(d, est_measure))
      do.call(
        rbind,
        lapply(probs, function(p) {
          lo <- (1 - p) / 2
          qs <- as.data.frame(posterior::summarise_draws(
            d,
            ~ posterior::quantile2(.x, probs = c(lo, 1 - lo))
          ))
          data.frame(
            row = i,
            category = qs$variable,
            estimate = est[[est_measure]],
            lower = qs[[2]],
            upper = qs[[3]],
            width = p
          )
        })
      )
    })
  )
  names(out)[names(out) == "category"] <- cat_col
  rownames(out) <- NULL

  # Ordered factor (widest first) so the default linewidth scale draws
  # narrower intervals more prominently; users can override via
  # `+ ggplot2::scale_linewidth_ordinal(...)`.
  width_labels <- sprintf("%g%%", probs * 100)
  out$width <- factor(sprintf("%g%%", out$width * 100), levels = width_labels)
  out
}

#' Build the `type = "interval"` plot for `plot_proportions()`: a
#' forest/caterpillar plot of one or more nested credible intervals per
#' observation, coloured by source and dodged along the x axis (rather
#' than faceted, so it scales to many observations without exploding into
#' tiny facets, and stays easy to further adjust with `ggplot2` calls),
#' narrower intervals drawn with a thicker `linewidth`. `group = source` is
#' set explicitly so both interval widths of the same source dodge
#' together, rather than `linewidth`'s levels (also a discrete aesthetic)
#' splitting each source into additional, separately-dodged groups.
#'
#' @param p_arr Numeric `[n_draws, n_obs, K]` array of proportion draws,
#'   with source names as the 3rd dimension's `dimnames`.
#' @param probs,robust,point_size See [plot_proportions()].
#' @param ... Further arguments passed to [ggplot2::geom_linerange()].
#' @return A `ggplot` object.
#' @noRd
plot_proportions_interval <- function(p_arr, probs, robust, point_size, ...) {
  df <- summarise_multi_interval(p_arr, probs, robust, cat_col = "source")
  df$row <- factor(df$row)
  dodge <- ggplot2::position_dodge(width = 0.5)

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$row, color = .data$source, group = .data$source)
  ) +
    ggplot2::geom_linerange(
      ggplot2::aes(
        ymin = .data$lower,
        ymax = .data$upper,
        linewidth = .data$width
      ),
      position = dodge,
      ...
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = .data$estimate),
      size = point_size,
      position = dodge
    ) +
    ggplot2::scale_linewidth_ordinal(range = c(0.5, 1.5)) +
    ggplot2::labs(
      x = "Observation",
      y = "Posterior proportion",
      color = "Source",
      linewidth = "Interval"
    ) +
    ggplot2::theme_minimal()
}
