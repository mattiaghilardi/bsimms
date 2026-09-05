#' Conditional effects of covariates on source proportions or predicted
#' mixture isotope values
#'
#' Computes posterior source proportions (`method = "posterior_proportions"`,
#' the default) or expected/predicted mixture isotope values (`method =
#' "posterior_epred"`/`"posterior_predict"`) across the range (for a numeric
#' covariate) or levels (for a factor covariate) of one or two
#' fixed-effect covariates from `formula`, holding every other
#' fixed-effect covariate at a reference value (the mean, for numeric
#' covariates; the first level, for factors -- overridable via
#' `ref_conditions`) and every group-level (random-effect) term at the
#' population-average level
#' (`re_formula = NA` by default, since the grid never varies grouping
#' columns). Printing or [plot()]-ing the result (requires `ggplot2`)
#' draws one plot per requested effect: a ribbon (numeric covariate) or
#' forest-plot-style linerange (factor covariate), with one or more
#' nested credible intervals (narrower drawn more prominently), faceted
#' by the second covariate when a two-way interaction is requested. For
#' `method = "posterior_proportions"`, every source is shown together on
#' one plot (colour/fill); for `"posterior_epred"`/`"posterior_predict"`,
#' isotopes are shown one at a time instead (see `resp`), since isotopes
#' don't share a common scale the way proportions do.
#'
#' @param object A `bsimms_fit` object (as returned by [bsimm()]).
#' @param x An object returned by `conditional_effects()`.
#' @param effects Character vector of fixed-effect covariate name(s) to
#'   vary (variable names in `formula`, not dummy-coded column names).
#'   `NULL` (default) computes conditional effects for every term of
#'   `formula`'s fixed-effect part, main effects and two-way interactions
#'   alike (e.g. `~ var1 * var2` gives conditional effects for `var1`,
#'   `var2`, and their interaction `var1:var2`); three-way and higher
#'   interactions are dropped from this default, since only two-way
#'   interactions are supported. An entry may name a
#'   two-way interaction as `"var1:var2"`: `var1` is swept as usual (its
#'   range if numeric, its observed levels if a factor) and becomes the
#'   plot's x axis, while `var2` becomes the *moderator* -- the plot is
#'   faceted over `var2`'s observed levels (factor) or a small set of
#'   representative values (numeric; see `int_conditions`).
#' @param ref_conditions Optional named list overriding the default
#'   reference value used to hold a fixed-effect covariate constant while
#'   another varies (any covariate currently neither swept nor used as a
#'   two-way interaction's moderator). Name each element after the
#'   covariate; the value is a single number (for a numeric covariate,
#'   replacing the default mean) or a single level (for a factor
#'   covariate, replacing the default of its first level). A covariate not
#'   named in `ref_conditions` uses the default.
#' @param int_conditions Optional named list overriding the default
#'   representative values used for a two-way interaction's moderator
#'   (`var2` in `"var1:var2"`; ignored for single covariates). Name each
#'   element after the moderator variable; the value is a numeric vector
#'   of representative values (for a numeric moderator, replacing the
#'   default mean and mean +/- 1 SD) or a character vector of levels to
#'   include (for a factor moderator, replacing the default of all
#'   observed levels). A moderator not named in `int_conditions` uses the
#'   default.
#' @param resolution Number of points spanning a numeric covariate's range
#'   (default 100). Ignored for factor covariates, which use their
#'   observed levels.
#' @param re_formula Which group-level terms to condition on; see
#'   [posterior_proportions()]. Defaults to `NA` (population-average),
#'   since the prediction grid only ever varies fixed-effect covariates.
#' @param robust Logical; if `FALSE` (default) the point estimate is the
#'   `mean`, if `TRUE` the `median`.
#' @param probs One or more credible-interval masses to display, e.g. the
#'   default `c(0.5, 0.95)` draws both a 50% and a 95% interval (nested,
#'   narrower intervals drawn more prominently).
#' @param point_size For factor covariates, the size of the point marking
#'   the central estimate (default 2). Ignored for numeric covariates.
#' @param method Which posterior quantity to compute and plot:
#'   `"posterior_proportions"` (default, via [posterior_proportions()]),
#'   `"posterior_epred"` (expected mixture isotope value, via
#'   [posterior_epred.bsimms_fit()]), or `"posterior_predict"` (posterior
#'   predictive mixture isotope value, via
#'   [posterior_predict.bsimms_fit()]); the latter two require
#'   `rstantools`.
#' @param resp For `method = "posterior_epred"`/`"posterior_predict"`,
#'   the isotope to plot (one of `isotope_names`); required if the model
#'   has more than one isotope, since isotopes are shown one at a time
#'   rather than together (unlike sources, they don't share a common
#'   scale). Defaults to the only isotope otherwise. Ignored for
#'   `method = "posterior_proportions"`.
#' @param plot For the `plot()` method, logical: display each plot as a
#'   side effect (default `TRUE`)? If `FALSE`, the plots are only built
#'   and returned.
#' @param ask For the `plot()` method, logical: prompt the user before
#'   displaying each new plot after the first (default `TRUE`)? Only
#'   relevant if `plot = TRUE` and there is more than one effect.
#' @param ... For `conditional_effects()`, further arguments passed to
#'   the underlying `method` function, e.g. `ndraws`. For the `plot()`
#'   method, further arguments passed to [ggplot2::geom_ribbon()]
#'   (numeric covariates) or [ggplot2::geom_linerange()] (factor
#'   covariates).
#' @return `conditional_effects()` returns an object of class
#'   `bsimms_conditional_effects`: a list of data frames, one per
#'   covariate/interaction in `effects`, each with one row per (grid
#'   point, category, interval width) and columns `row`, `source`/
#'   `isotope` (depending on `method`), `estimate`, `lower`, `upper`,
#'   `width`, the covariate's values, and (for a two-way interaction) the
#'   moderator's values. The `plot()` method returns a list of `ggplot`
#'   objects (one per covariate), invisibly; each is also displayed as a
#'   side effect.
#' @export
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
#' fit <- bsimm(
#'   sim$formula, mixture_data = sim$mixture_data,
#'   source_data = sim$source_data, tdf_data = sim$tdf_data,
#'   isotope_names = sim$isotope_names,
#'   source_means_sds = sim$source_means_sds, tdf_means_sds = sim$tdf_means_sds,
#'   conc_dep = sim$conc_dep, error_structure = sim$error_structure,
#'   source_col = sim$source_col,
#'   chains = 2, iter_warmup = 500, iter_sampling = 500
#' )
#' ce <- conditional_effects(fit)
#' plot(ce)
#' }
conditional_effects <- function(
  object,
  effects = NULL,
  ref_conditions = NULL,
  int_conditions = NULL,
  resolution = 100,
  re_formula = NA,
  robust = FALSE,
  probs = c(0.5, 0.95),
  point_size = 2,
  method = c("posterior_proportions", "posterior_epred", "posterior_predict"),
  resp = NULL,
  ...
) {
  if (!inherits(object, "bsimms_fit")) {
    cli::cli_abort(
      "{.arg object} must be a {.cls bsimms_fit} object.",
      call = NULL
    )
  }
  method <- rlang::arg_match(method)
  spec <- object$spec
  available <- names(spec$fixed_frame)
  if (length(available) == 0) {
    cli::cli_abort(
      "Model has no fixed-effect covariates to condition on.",
      call = NULL
    )
  }
  check_condition_names(ref_conditions, "ref_conditions", available)
  check_condition_names(int_conditions, "int_conditions", available)

  if (method != "posterior_proportions") {
    rlang::check_installed(
      "rstantools",
      reason = paste0(
        "to use `conditional_effects()` with method = ",
        "\"posterior_epred\"/\"posterior_predict\"."
      )
    )
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
  }
  cat_col <- if (method == "posterior_proportions") "source" else "isotope"
  y_label <- if (method == "posterior_proportions") {
    "Posterior proportion"
  } else {
    resp
  }

  effects <- if (is.null(effects)) {
    default_conditional_effects(spec$fixed_formula)
  } else {
    effects
  }
  effect_vars <- lapply(
    effects,
    parse_conditional_effect,
    available = available
  )
  names(effect_vars) <- effects

  out <- lapply(effect_vars, function(vars) {
    cond <- build_conditional_grid(
      spec$fixed_frame,
      vars,
      ref_conditions,
      int_conditions,
      resolution
    )
    arr <- switch(
      method,
      posterior_proportions = posterior_proportions(
        object,
        newdata = cond$grid,
        re_formula = re_formula,
        ...
      ),
      posterior_epred = rstantools::posterior_epred(
        object,
        newdata = cond$grid,
        resp = resp,
        re_formula = re_formula,
        ...
      ),
      posterior_predict = rstantools::posterior_predict(
        object,
        newdata = cond$grid,
        resp = resp,
        re_formula = re_formula,
        ...
      )
    )
    df <- summarise_multi_interval(arr, probs, robust, cat_col = cat_col)
    df[[vars[1]]] <- cond$grid[[vars[1]]][df$row]
    attr(df, "cat_col") <- cat_col
    attr(df, "y_label") <- y_label
    attr(df, "effect") <- vars[1]
    attr(df, "is_numeric") <- cond$is_numeric
    attr(df, "point_size") <- point_size
    if (length(vars) == 2) {
      df[[vars[2]]] <- cond$grid[[vars[2]]][df$row]
      attr(df, "moderator") <- vars[2]
      attr(df, "moderator_is_numeric") <- cond$moderator_is_numeric
    }
    df
  })
  class(out) <- "bsimms_conditional_effects"
  out
}

#' @rdname conditional_effects
#' @export
plot.bsimms_conditional_effects <- function(x, plot = TRUE, ask = TRUE, ...) {
  rlang::check_installed(
    "ggplot2",
    reason = "to use `plot.bsimms_conditional_effects()`."
  )
  plots <- lapply(x, plot_one_conditional_effect, ...)
  names(plots) <- names(x)
  if (plot) {
    default_ask <- grDevices::devAskNewPage()
    on.exit(grDevices::devAskNewPage(default_ask))
    grDevices::devAskNewPage(ask = FALSE)
    for (i in seq_along(plots)) {
      print(plots[[i]])
      if (i == 1) grDevices::devAskNewPage(ask = ask)
    }
  }
  invisible(plots)
}

#' @export
print.bsimms_conditional_effects <- function(x, ...) {
  plot(x, ...)
}

#' Default `effects` (when the user supplies `NULL`): every term of
#' `fixed_formula`, main effects and two-way interactions alike, e.g.
#' `~ var1 * var2` gives `c("var1", "var2", "var1:var2")`. Three-way and
#' higher interactions are dropped rather than included (and erroring,
#' since `conditional_effects()` only supports two-way interactions) -- a
#' model with one still gets conditional effects for its lower-order
#' terms by default; the higher-order term itself must be requested some
#' other way (not currently supported).
#'
#' @param fixed_formula The model's fixed-effect formula
#'   (`spec$fixed_formula`).
#' @return Character vector of term labels (see `stats::terms()`), with
#'   any 3+-way interaction terms removed.
#' @noRd
default_conditional_effects <- function(fixed_formula) {
  term_labels <- attr(stats::terms(fixed_formula), "term.labels")
  Filter(function(term) length(strsplit(term, ":")[[1]]) <= 2, term_labels)
}

#' Validate that every name in `int_conditions`/`ref_conditions` refers to
#' an actual fixed-effect covariate in the model, catching typos that
#' would otherwise be silently ignored: `x[[name]]` for a `name` that
#' isn't really one of `x`'s keys (e.g. a misspelled covariate) simply
#' returns `NULL`, so the lookup falls back to the default instead of
#' erroring.
#'
#' @param x `int_conditions` or `ref_conditions` (see
#'   [conditional_effects()]), or `NULL`.
#' @param arg_name `"int_conditions"` or `"ref_conditions"`, for the error
#'   message.
#' @param available Character vector of the model's fixed-effect covariate
#'   names (`names(spec$fixed_frame)`).
#' @return Invisible `NULL`; called for its error side effect.
#' @noRd
check_condition_names <- function(x, arg_name, available) {
  unknown <- setdiff(names(x), available)
  if (length(unknown) > 0) {
    cli::cli_abort(
      paste0(
        "{.arg {arg_name}} names covariate{?s} not in the model: ",
        "{.val {unknown}}."
      ),
      call = NULL
    )
  }
}

#' Parse one `effects` entry into one or two validated covariate names.
#'
#' @param eff A single `effects` entry: a covariate name, or `"var1:var2"`
#'   naming a two-way interaction.
#' @param available Character vector of the model's fixed-effect covariate
#'   names (`names(spec$fixed_frame)`).
#' @return A character vector of length 1 or 2: the validated covariate
#'   name(s) in `eff`.
#' @noRd
parse_conditional_effect <- function(eff, available) {
  parts <- trimws(strsplit(eff, ":")[[1]])
  if (length(parts) > 2) {
    cli::cli_abort(
      paste0(
        "{.arg effects} entry {.val {eff}} names more than two ",
        "covariates; only two-way interactions are supported."
      ),
      call = NULL
    )
  }
  unknown <- setdiff(parts, available)
  if (length(unknown) > 0) {
    cli::cli_abort(
      paste0(
        "{.arg effects} entry {.val {eff}} refers to covariate{?s} not ",
        "in the model: {.val {unknown}}."
      ),
      call = NULL
    )
  }
  if (length(parts) == 2 && parts[1] == parts[2]) {
    cli::cli_abort(
      "{.arg effects} entry {.val {eff}} names the same covariate twice.",
      call = NULL
    )
  }
  parts
}

#' Resolve a two-way interaction's moderator (`vars[2]`) representative
#' values: the matching element of `int_conditions` if supplied (validated
#' -- numeric for a numeric moderator, a subset of its levels for a
#' factor moderator), otherwise the default (all observed levels for a
#' factor moderator, mean and mean +/- 1 SD for a numeric one).
#'
#' @param mod_col The moderator's column in `fixed_frame`.
#' @param moderator The moderator's name (for error messages and looking
#'   it up in `int_conditions`).
#' @param int_conditions See [conditional_effects()].
#' @return A numeric or character vector of representative values.
#' @noRd
resolve_moderator_values <- function(mod_col, moderator, int_conditions) {
  is_numeric <- is.numeric(mod_col)
  values <- int_conditions[[moderator]]

  if (is.null(values)) {
    return(
      if (is_numeric) {
        mean(mod_col, na.rm = TRUE) +
          c(-1, 0, 1) * stats::sd(mod_col, na.rm = TRUE)
      } else {
        levels(mod_col)
      }
    )
  }
  if (is_numeric) {
    if (!is.numeric(values)) {
      cli::cli_abort(
        paste0(
          "{.arg int_conditions} for {.field {moderator}} must be ",
          "numeric, since it is a numeric covariate."
        ),
        call = NULL
      )
    }
    return(values)
  }
  values <- as.character(values)
  unknown <- setdiff(values, levels(mod_col))
  if (length(unknown) > 0) {
    cli::cli_abort(
      c(
        paste0(
          "{.arg int_conditions} for {.field {moderator}} contains ",
          "level{?s} not in the model: {.val {unknown}}."
        ),
        "i" = paste0(
          "Valid levels for {.field {moderator}} are: ",
          "{.val {levels(mod_col)}}."
        )
      ),
      call = NULL
    )
  }
  values
}

#' Resolve a fixed-effect covariate's reference value (used to hold it
#' constant while another covariate varies): the matching element of
#' `ref_conditions` if supplied (validated -- a single number for a
#' numeric covariate, a single valid level for a factor covariate),
#' otherwise the default (the mean for a numeric covariate, its first
#' level for a factor one).
#'
#' @param col The covariate's column in `fixed_frame`.
#' @param name The covariate's name (for error messages and looking it up
#'   in `ref_conditions`).
#' @param ref_conditions See [conditional_effects()].
#' @return A single number or factor value.
#' @noRd
resolve_reference_value <- function(col, name, ref_conditions) {
  is_numeric <- is.numeric(col)
  value <- ref_conditions[[name]]

  if (is.null(value)) {
    return(
      if (is_numeric) {
        mean(col, na.rm = TRUE)
      } else {
        factor(levels(col)[1], levels = levels(col))
      }
    )
  }
  if (length(value) != 1) {
    cli::cli_abort(
      "{.arg ref_conditions} for {.field {name}} must be a single value.",
      call = NULL
    )
  }
  if (is_numeric) {
    if (!is.numeric(value)) {
      cli::cli_abort(
        paste0(
          "{.arg ref_conditions} for {.field {name}} must be numeric, ",
          "since it is a numeric covariate."
        ),
        call = NULL
      )
    }
    return(value)
  }
  value <- as.character(value)
  if (!value %in% levels(col)) {
    cli::cli_abort(
      c(
        paste0(
          "{.arg ref_conditions} for {.field {name}} contains a level ",
          "not in the model: {.val {value}}."
        ),
        "i" = "Valid levels for {.field {name}} are: {.val {levels(col)}}."
      ),
      call = NULL
    )
  }
  factor(value, levels = levels(col))
}

#' Build the plot for one covariate's conditional effect: a ribbon (numeric
#' covariate) or forest-plot-style linerange (factor covariate), with one
#' or more nested credible intervals (narrower drawn more prominently --
#' `alpha` for ribbons, `linewidth` for lineranges), faceted by the
#' moderator when `df` carries one (see `conditional_effects()`'s two-way
#' interaction support). Coloured by category (`cat_col`, `"source"` or
#' `"isotope"`) when `df` has more than one -- true for
#' `method = "posterior_proportions"` (every source shown together), but
#' not for `"posterior_epred"`/`"posterior_predict"` (exactly one isotope,
#' via `resp`), which are plotted plain, with no legend.
#'
#' @param df One element of a `bsimms_conditional_effects` object (see
#'   [conditional_effects()]): a data frame with columns `row`,
#'   `<cat_col>`, `estimate`, `lower`, `upper`, `width`, the covariate's
#'   values, and (for a two-way interaction) the moderator's values, plus
#'   attributes `effect` (the covariate name), `is_numeric`, `point_size`,
#'   `cat_col`, `y_label`, and (if a two-way interaction)
#'   `moderator`/`moderator_is_numeric`.
#' @param ... Further arguments passed to [ggplot2::geom_ribbon()] or
#'   [ggplot2::geom_linerange()]; see [plot.bsimms_conditional_effects()].
#' @return A `ggplot` object.
#' @noRd
plot_one_conditional_effect <- function(df, ...) {
  effect <- attr(df, "effect")
  is_numeric <- attr(df, "is_numeric")
  point_size <- attr(df, "point_size")
  moderator <- attr(df, "moderator")
  cat_col <- attr(df, "cat_col")
  y_label <- attr(df, "y_label")
  show_legend <- length(unique(df[[cat_col]])) > 1
  legend_label <- if (cat_col == "source") "Source" else "Isotope"

  facet_layer <- NULL
  if (!is.null(moderator)) {
    if (attr(df, "moderator_is_numeric")) {
      mod_label <- paste0(moderator, ": ", signif(df[[moderator]], 3))
    } else {
      mod_label <- as.character(df[[moderator]])
    }
    df$.moderator_label <- factor(
      mod_label,
      levels = unique(mod_label[order(df[[moderator]])])
    )
    facet_layer <- ggplot2::facet_wrap(~.moderator_label)
  }

  if (is_numeric) {
    line_cols <- c(
      effect,
      cat_col,
      "estimate",
      if (!is.null(moderator)) ".moderator_label"
    )
    line_df <- unique(df[, line_cols])
    mapping <- if (show_legend) {
      ggplot2::aes(
        x = .data[[effect]],
        color = .data[[cat_col]],
        fill = .data[[cat_col]]
      )
    } else {
      ggplot2::aes(x = .data[[effect]])
    }
    g <- ggplot2::ggplot(df, mapping) +
      ggplot2::geom_ribbon(
        ggplot2::aes(
          ymin = .data$lower,
          ymax = .data$upper,
          alpha = .data$width
        ),
        color = NA,
        ...
      ) +
      ggplot2::geom_line(data = line_df, ggplot2::aes(y = .data$estimate)) +
      ggplot2::scale_alpha_ordinal(range = c(0.4, 0.15)) +
      facet_layer +
      ggplot2::labs(x = effect, y = y_label, alpha = "Interval") +
      ggplot2::theme_minimal()
    if (show_legend) {
      g <- g + ggplot2::labs(color = legend_label, fill = legend_label)
    }
    g
  } else {
    dodge <- ggplot2::position_dodge(width = 0.5)
    mapping <- if (show_legend) {
      ggplot2::aes(
        x = .data[[effect]],
        color = .data[[cat_col]],
        group = .data[[cat_col]]
      )
    } else {
      # `group` must still be set, even with a single category: without it,
      # `position_dodge()` falls back to grouping by any discrete aesthetic
      # in scope -- including `linewidth = width` (mapped below) -- which
      # would dodge the nested interval widths apart instead of overlaying
      # them at the same x.
      ggplot2::aes(x = .data[[effect]], group = .data[[cat_col]])
    }
    g <- ggplot2::ggplot(df, mapping) +
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
      facet_layer +
      ggplot2::labs(x = effect, y = y_label, linewidth = "Interval") +
      ggplot2::theme_minimal()
    if (show_legend) {
      g <- g + ggplot2::labs(color = legend_label)
    }
    g
  }
}

#' Build a prediction grid for [conditional_effects()]: `vars[1]` spans its
#' range (numeric) or observed levels (factor) and becomes the plot's x
#' axis; if `vars` also names a moderator (`vars[2]`), the grid crosses
#' `vars[1]`'s sweep with the moderator's representative values (see
#' `resolve_moderator_values()`) via `expand.grid()`. Every other
#' fixed-effect covariate in `fixed_frame` is held at a reference value
#' (see `resolve_reference_value()`).
#'
#' @param fixed_frame Model frame of the fixed-effect covariates (as built
#'   by `parse_bsimms_formula()`, i.e. `spec$fixed_frame`).
#' @param vars Character vector of length 1 or 2 (as returned by
#'   `parse_conditional_effect()`): the covariate(s) to vary.
#' @param ref_conditions,int_conditions See [conditional_effects()].
#' @param resolution Number of points spanning `vars[1]`'s range, if
#'   numeric. Ignored for a factor `vars[1]`, which uses its observed
#'   levels.
#' @return A list with `grid` (a data frame, one row per grid point, same
#'   columns as `fixed_frame`), `is_numeric` (logical, whether `vars[1]`
#'   is a numeric covariate), and `moderator_is_numeric` (logical, or
#'   `NULL` if `vars` has no moderator).
#' @noRd
build_conditional_grid <- function(
  fixed_frame,
  vars,
  ref_conditions,
  int_conditions,
  resolution
) {
  effect <- vars[1]
  focal <- fixed_frame[[effect]]
  is_numeric <- is.numeric(focal)
  focal_values <- if (is_numeric) {
    seq(
      min(focal, na.rm = TRUE),
      max(focal, na.rm = TRUE),
      length.out = resolution
    )
  } else {
    levels(focal)
  }

  moderator_is_numeric <- NULL
  if (length(vars) == 2) {
    moderator <- vars[2]
    mod_col <- fixed_frame[[moderator]]
    moderator_is_numeric <- is.numeric(mod_col)
    moderator_values <- resolve_moderator_values(
      mod_col,
      moderator,
      int_conditions
    )
    grid <- expand.grid(
      focal_values,
      moderator_values,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    names(grid) <- vars
  } else {
    grid <- data.frame(x = focal_values)
    names(grid) <- effect
  }

  for (v in setdiff(names(fixed_frame), vars)) {
    ref <- resolve_reference_value(fixed_frame[[v]], v, ref_conditions)
    grid[[v]] <- rep(ref, nrow(grid))
  }
  list(
    grid = grid,
    is_numeric = is_numeric,
    moderator_is_numeric = moderator_is_numeric
  )
}
