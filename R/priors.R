#' Specify a prior for one or more `bsimms` parameters
#'
#' Build up a full prior specification by combining several calls with
#' [c()], e.g.
#' `c(bsimms_prior("normal(0, 2)", class = "b"),
#'    bsimms_prior("normal(1, 1)", class = "b", coef = "SeasonWinter"))`.
#' More specific rows (a given `coef`/`resp`/`group`) take precedence over
#' the general class default when the Stan code is generated.
#'
#' @param prior Character string: a valid Stan distribution expression,
#'   e.g. `"normal(0, 1)"`, `"student_t(3, 0, 2.5)"`, `"lkj_corr_cholesky(2)"`,
#'   or, for `class = "p_global"`, a single positive number (a Dirichlet
#'   concentration, e.g. `"1"`).
#' @param class One of:
#'   - `"b"`: fixed-effect slopes (the population-level baseline is not a
#'     `"b"` coefficient but `"p_global"`, see below).
#'   - `"p_global"`: Dirichlet concentration for one source's share of the
#'     global/population-average proportions (MixSIAR's `p_global`).
#'   - `"sd"`: group-level standard deviations.
#'   - `"cor"`: group-level correlations (LKJ prior on the Cholesky
#'     factor).
#'   - `"sigma"`: residual/observation error; only used for
#'     `error_structure` `"residual_only"`.
#'   - `"resid_prop"`: MixSIAR's `resid.prop`, a multiplicative factor
#'     scaling the propagated source/TDF process variance; only used for
#'     `error_structure` `"process_residual"`; always restricted to
#'     values between 0 and 20, matching the range of the recommended
#'     default prior (see [bsimms_get_prior()]), so a custom prior can
#'     reshape which values in that range are more likely but cannot
#'     allow values above 20.
#'   - `"source_mean"`, `"source_sd"`, `"tdf_mean"`, `"tdf_sd"`: only
#'     used when the corresponding data are supplied raw rather than as
#'     means/SDs.
#'   - `"source_cor"`: LKJ prior on the Cholesky factor of each source's
#'     isotope correlation matrix; only used when source data are raw
#'     and there are 2+ isotopes.
#'   - `"resid_cor"`: LKJ prior on the Cholesky factor of the shared
#'     residual-error correlation matrix; only used for
#'     `error_structure` `"residual_only"` with 2+ isotopes.
#' @param coef Optional: restrict to one fixed-effect coefficient name, as
#'   it appears in the model formula's expanded design matrix (e.g.
#'   `"SexM"` for a factor `Sex`, or `"SexM:RegionB"` for an interaction).
#'   Only used with `class = "b"`.
#' @param resp Optional: restrict to one isotope (response) name. Used with
#'   `class` `"sigma"`, `"resid_prop"`, `"source_mean"`, `"source_sd"`,
#'   `"tdf_mean"` or `"tdf_sd"`.
#' @param group Optional, one of two uses depending on `class`:
#'   - For `"sd"`/`"cor"`, restrict to one group-level term (the
#'     right-hand side of a `(... | group)` term, as written in the
#'     formula).
#'   - For `"p_global"`/`"source_mean"`/`"source_sd"`/`"tdf_mean"`/
#'     `"tdf_sd"`/`"source_cor"`, restrict to one source (as named in
#'     `source_data`/`tdf_data`), so that source and trophic
#'     discrimination factor priors can be made source-specific as well
#'     as isotope-specific. Left unset (`""`), such a prior applies to
#'     every source (and, except for `"source_cor"`, every isotope).
#' @return A one-row `data.frame` of class `bsimms_prior`.
#' @export
#' @examples
#' bsimms_prior("normal(0, 2)", class = "b")
#' bsimms_prior("normal(1, 1)", class = "b", coef = "SeasonWinter")
#' bsimms_prior("student_t(3, 0, 1)", class = "sd", group = "Region")
#' bsimms_prior(
#'   "normal(3.4, 0.3)",
#'   class = "tdf_mean", resp = "d15N", group = "Beaver"
#' )
#' bsimms_prior(
#'   "student_t(3, 0, 0.5)",
#'   class = "source_sd", resp = "d13C", group = "Beaver"
#' )
#' bsimms_prior("2", class = "p_global", group = "Beaver")
#' bsimms_prior("lkj_corr_cholesky(2)", class = "source_cor", group = "Beaver")
#' bsimms_prior("lkj_corr_cholesky(2)", class = "resid_cor")
bsimms_prior <- function(prior, class = "b", coef = "", resp = "", group = "") {
  class <- rlang::arg_match0(class, c(
    "b", "p_global", "sd", "cor", "sigma", "resid_prop",
    "source_mean", "source_sd", "tdf_mean", "tdf_sd", "source_cor", "resid_cor"
  ))
  if (!rlang::is_string(prior)) {
    cli::cli_abort(
      "{.arg prior} must be a single character string of Stan code, e.g. {.code \"normal(0, 1)\"}.",
      call = NULL
    )
  }
  out <- data.frame(
    prior = prior, class = class, coef = coef, resp = resp, group = group,
    stringsAsFactors = FALSE
  )
  class(out) <- c("bsimms_prior", "data.frame")
  out
}

#' @export
c.bsimms_prior <- function(...) {
  out <- do.call(rbind, lapply(list(...), as.data.frame))
  class(out) <- c("bsimms_prior", "data.frame")
  out
}

#' Print a `bsimms_prior` specification
#'
#' Prints as a table with columns `prior`, `class`, `coef`, `resp` and
#' `group` (blank where unset).
#'
#' @param x A `bsimms_prior` object (as returned by [bsimms_prior()] or
#'   [bsimms_get_prior()]).
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.bsimms_prior <- function(x, ...) {
  print.data.frame(as.data.frame(x), row.names = FALSE, quote = FALSE)
  invisible(x)
}

#' Build the table of default priors for a model, as returned by
#' [bsimms_get_prior()]: one class-level row per parameter class the model
#' actually has (given `spec`'s dimensions, `error_structure`, and whether
#' source/TDF data are raw or summarised), plus one `group`-specific row per
#' source (`"p_global"`), one `resp`-specific row per isotope (`"sigma"`/
#' `"resid_prop"`, where applicable), and one `resp`- and `group`-specific
#' row per isotope/source combination (`"source_mean"`/`"source_sd"`/
#' `"tdf_mean"`/`"tdf_sd"`, since these are estimated separately per source
#' when the corresponding data are raw), each pre-filled with a weakly
#' informative default. Where a natural data scale exists (`sigma`/
#' `source_mean`/`source_sd`/`tdf_mean`/`tdf_sd`), the default is scaled
#' using the sample median/MAD (median absolute deviation) of the relevant
#' data -- for `source_mean`/`source_sd`/`tdf_mean`/`tdf_sd`, that source's
#' own raw replicates only, so the prior reflects that source's actual
#' location/dispersion rather than being diluted by other sources: median
#' and MAD are used instead of the mean and SD because they are robust to
#' outliers and skew, which isotope data are prone to.
#'
#' @param spec A `bsimms_spec` (see `build_bsimms_spec()`).
#' @return A `bsimms_prior` data frame (see [bsimms_prior()]).
#' @noRd
default_bsimms_prior <- function(spec) {
  rows <- list()
  add <- function(prior, class, coef = "", resp = "", group = "") {
    rows[[length(rows) + 1]] <<- bsimms_prior(prior, class, coef, resp, group)
  }

  add("normal(0, 1)", "b")

  for (s in spec$source_names) add("1", "p_global", group = s)

  if (length(spec$re_terms) > 0) {
    add("student_t(3, 0, 1)", "sd")
    add("lkj_corr_cholesky(1)", "cor")
  }

  if (spec$error_structure == "residual_only") {
    y_mad <- apply(spec$y, 2, stats::mad)
    y_mad[!is.finite(y_mad) | y_mad <= 0] <- 1
    add(sprintf("student_t(3, 0, %.6g)", max(y_mad)), "sigma")
    for (j in seq_along(spec$isotope_names)) {
      add(sprintf("student_t(3, 0, %.6g)", y_mad[j]), "sigma", resp = spec$isotope_names[j])
    }
    if (spec$J > 1) {
      add("lkj_corr_cholesky(1)", "resid_cor")
    }
  }

  if (spec$error_structure == "process_residual") {
    # uniform(0, 20), not a data-scaled default: Stock & Semmens (2016)
    # Appendix S2 simulation-tested it against three priors centred/moded at
    # the expected value 1 (chi-square(3), gamma(2,2), log-normal(0,1)) and
    # found it better-calibrated at high consumer variance, so it is kept
    # exactly as recommended rather than replaced with a more "natural"-
    # looking peaked prior.
    add("uniform(0, 20)", "resid_prop")
    for (j in seq_along(spec$isotope_names)) {
      add("uniform(0, 20)", "resid_prop", resp = spec$isotope_names[j])
    }
  }

  if (spec$source$mode == "raw") {
    Y <- spec$source$Y
    idx <- spec$source$source_idx
    for (j in seq_along(spec$isotope_names)) {
      for (k in seq_along(spec$source_names)) {
        yk <- Y[idx == k, j]
        m <- stats::median(yk); s <- max(stats::mad(yk), 1e-3)
        add(
          sprintf("normal(%.6g, %.6g)", m, 10 * s), "source_mean",
          resp = spec$isotope_names[j], group = spec$source_names[k]
        )
        add(
          sprintf("student_t(3, 0, %.6g)", s), "source_sd",
          resp = spec$isotope_names[j], group = spec$source_names[k]
        )
      }
    }
    if (spec$J > 1) {
      add("lkj_corr_cholesky(1)", "source_cor")
      for (k in seq_along(spec$source_names)) {
        add("lkj_corr_cholesky(1)", "source_cor", group = spec$source_names[k])
      }
    }
  }
  if (spec$tdf$mode == "raw") {
    Y <- spec$tdf$Y
    idx <- spec$tdf$source_idx
    for (j in seq_along(spec$isotope_names)) {
      for (k in seq_along(spec$source_names)) {
        yk <- Y[idx == k, j]
        m <- stats::median(yk); s <- max(stats::mad(yk), 1e-3)
        add(
          sprintf("normal(%.6g, %.6g)", m, 10 * s), "tdf_mean",
          resp = spec$isotope_names[j], group = spec$source_names[k]
        )
        add(
          sprintf("student_t(3, 0, %.6g)", s), "tdf_sd",
          resp = spec$isotope_names[j], group = spec$source_names[k]
        )
      }
    }
  }

  do.call(c, rows)
}

#' Validate that a user-supplied prior's `coef`/`resp`/`group` values
#' (where meaningful for that row's `class`) are real identifiers from
#' `spec`, rather than a typo that would otherwise be silently appended as
#' a dead row -- never matched by `select_prior()`, so the override would
#' never actually apply to anything and no error would ever surface it.
#'
#' @param user A `bsimms_prior` data frame of user overrides.
#' @param spec A `bsimms_spec` (see `build_bsimms_spec()`).
#' @return Invisible `NULL`; errors on an unrecognised identifier.
#' @noRd
validate_prior_identifiers <- function(user, spec) {
  re_groups <- vapply(spec$re_terms, function(x) x$group, character(1))
  check <- function(value, valid, field, class) {
    if (!nzchar(value) || value %in% valid) {
      return(invisible(NULL))
    }
    if (length(valid) == 0) {
      cli::cli_abort(
        "{.arg prior}: class {.val {class}} has no valid {.field {field}} in this model (no matching terms).",
        call = NULL
      )
    }
    cli::cli_abort(
      "{.arg prior}: unrecognised {.field {field}} {.val {value}} for class {.val {class}}; must be one of {.val {valid}}.",
      call = NULL
    )
  }
  for (i in seq_len(nrow(user))) {
    class <- user$class[i]
    if (class == "b") check(user$coef[i], spec$fixed_names, "coef", class)
    if (class %in% c("sigma", "resid_prop", "source_mean", "source_sd", "tdf_mean", "tdf_sd")) {
      check(user$resp[i], spec$isotope_names, "resp", class)
    }
    if (class %in% c("sd", "cor")) check(user$group[i], re_groups, "group", class)
    if (class %in% c("p_global", "source_mean", "source_sd", "tdf_mean", "tdf_sd", "source_cor")) {
      check(user$group[i], spec$source_names, "group", class)
    }
  }
  invisible(NULL)
}

#' Merge a user-supplied prior specification into the model's default
#' priors: rows matching an existing `(class, coef, resp, group)`
#' combination overwrite that row's `prior` string in place, and rows with
#' no match are appended (allowing new, more specific overrides, e.g. a
#' `resp`-specific row when only a class-level default exists). Finally,
#' any user row that is *less* specific than existing rows of the same
#' class (e.g. a class-level `"resid_prop"` override, when
#' `default_bsimms_prior()` also pre-fills a `resp`-specific row for every
#' isotope) is cascaded onto those more specific rows too -- but never
#' onto one the user separately overrode themselves, which keeps its own,
#' more specific value regardless of the two rows' order. Without this,
#' `select_prior()`'s always-most-specific-first lookup would find the
#' untouched, pre-filled specific row and the less specific override
#' would silently never be used.
#'
#' @param default A `bsimms_prior` data frame (as returned by
#'   `default_bsimms_prior()`).
#' @param user `NULL`, or a `bsimms_prior` data frame of user overrides (as
#'   returned by [bsimms_prior()], optionally combined with [c()]).
#' @param spec A `bsimms_spec`, used to validate `user`'s `coef`/`resp`/
#'   `group` values against the model's real identifiers (see
#'   `validate_prior_identifiers()`).
#' @return A `bsimms_prior` data frame: `default` unchanged if `user` is
#'   `NULL`, otherwise `default` with `user`'s rows merged in.
#' @noRd
merge_bsimms_prior <- function(default, user, spec) {
  if (is.null(user)) return(default)
  if (!inherits(user, "bsimms_prior")) {
    cli::cli_abort(
      "{.arg prior} must be built with {.fn bsimms_prior} (optionally combined with {.fn c}).",
      call = NULL
    )
  }
  validate_prior_identifiers(user, spec)
  out <- default
  touched <- rep(FALSE, nrow(out))
  for (i in seq_len(nrow(user))) {
    match_row <- which(
      out$class == user$class[i] & out$coef == user$coef[i] &
        out$resp == user$resp[i] & out$group == user$group[i]
    )
    if (length(match_row) == 1) {
      out$prior[match_row] <- user$prior[i]
      touched[match_row] <- TRUE
    } else {
      out <- rbind(out, user[i, ])
      touched <- c(touched, TRUE)
    }
  }
  for (i in seq_len(nrow(user))) {
    candidates <- which(out$class == user$class[i] & !touched)
    for (r in candidates) {
      if (is_more_specific_prior_row(out[r, ], user[i, ])) {
        out$prior[r] <- user$prior[i]
      }
    }
  }
  class(out) <- c("bsimms_prior", "data.frame")
  out
}

#' Is `row` a specialisation of `base` -- i.e. do they agree on every
#' `coef`/`resp`/`group` field where `base` is non-empty, and does `row`
#' additionally pin down at least one field `base` leaves unrestricted
#' (`""`)? Used by `merge_bsimms_prior()` to cascade a less specific user
#' override onto the more specific default rows it would otherwise be
#' shadowed by, mirroring `select_prior()`'s own specificity ordering.
#'
#' @param row,base Single-row `bsimms_prior` data frames (same `class`).
#' @return Logical.
#' @noRd
is_more_specific_prior_row <- function(row, base) {
  fields <- c("coef", "resp", "group")
  more_specific <- FALSE
  for (f in fields) {
    if (nzchar(base[[f]])) {
      if (row[[f]] != base[[f]]) return(FALSE)
    } else if (nzchar(row[[f]])) {
      more_specific <- TRUE
    }
  }
  more_specific
}

#' Look up the Stan prior expression to use for one parameter, from a
#' (merged) prior table. Falls back from most to least specific: tries
#' every combination of `(coef, resp, group)` as given down to the fully
#' unrestricted class-level default, in decreasing order of how many of
#' the three are used; ties within a specificity level are broken by
#' taking the last matching row (so a later `merge_bsimms_prior()`
#' override always wins over an earlier one at the same specificity).
#'
#' @param prior_df A `bsimms_prior` data frame (typically the merged
#'   default + user table).
#' @param class One of `bsimms_prior()`'s `class` values.
#' @param coef,resp,group Optional `coef`/`resp`/`group` to match, as in
#'   [bsimms_prior()]. `""` (default) matches the class-level default.
#' @return A single character string: the matching Stan prior expression.
#' @noRd
select_prior <- function(prior_df, class, coef = "", resp = "", group = "") {
  sub <- prior_df[prior_df$class == class, , drop = FALSE]
  if (nrow(sub) == 0) {
    cli::cli_abort(
      "No prior available for class {.val {class}}; this should not happen \u2014 please report a bug.",
      call = NULL
    )
  }
  pick <- function(co, rs, gr) {
    hit <- which(sub$coef == co & sub$resp == rs & sub$group == gr)
    if (length(hit) >= 1) utils::tail(sub$prior[hit], 1) else NA_character_
  }
  candidates <- list(
    pick(coef, resp, group),
    pick(coef, resp, ""), pick(coef, "", group), pick("", resp, group),
    pick(coef, "", ""), pick("", resp, ""), pick("", "", group),
    pick("", "", "")
  )
  candidates <- candidates[!vapply(candidates, is.na, logical(1))]
  if (length(candidates) == 0) {
    cli::cli_abort(
      "No matching prior found for class {.val {class}}, coef {.val {coef}}, resp {.val {resp}}, group {.val {group}}.",
      call = NULL
    )
  }
  candidates[[1]]
}

#' Default priors for a `bsimms` model
#'
#' Builds the model specification (design matrices, source/TDF data layout)
#' without generating or fitting the Stan model, and returns the table of
#' default priors that would be used. Edit
#' rows of the result and pass to the `prior` argument of `make_stancode()`,
#' `make_standata()` or `bsimm()` to override.
#'
#' @inheritParams make_stancode
#' @return A `bsimms_prior` data frame.
#' @references Moore, J.W., & Semmens, B.X. (2008). Incorporating
#'   uncertainty and prior information into stable isotope mixing models.
#'   *Ecology Letters*, 11(5), 470-480.
#'   \doi{10.1111/j.1461-0248.2008.01163.x}
#' @references Parnell, A.C., Inger, R., Bearhop, S., & Jackson, A.L.
#'   (2010). Source partitioning using stable isotopes: coping with too
#'   much variation. *PLoS ONE*, 5(3), e9672.
#'   \doi{10.1371/journal.pone.0009672}
#' @references Stock, B.C., & Semmens, B.X. (2016). Unifying error
#'   structures in commonly used biotracer mixing models. *Ecology*,
#'   97(10), 2562-2569. \doi{10.1002/ecy.1517} Appendix S2 simulation-tests
#'   the default `"resid_prop"` prior (`uniform(0, 20)`) against three
#'   priors centred/moded at 1 (chi-square(3), gamma(2,2), log-normal(0,1))
#'   and finds it better-calibrated at high consumer variance.
#' @export
#' @examples
#' sim <- simulate_bsimms_data(
#'   ~ 1 + (1 | Region),
#'   n_mixture_obs = 20,
#'   source_names = c("Beaver", "Deer", "Hare"),
#'   isotope_names = c("d13C", "d15N"),
#'   n_groups = list(Region = 2),
#'   seed = 1
#' )
#' bsimms_get_prior(
#'   sim$formula, mixture_data = sim$mixture_data,
#'   source_data = sim$source_data, tdf_data = sim$tdf_data,
#'   isotope_names = sim$isotope_names,
#'   source_means_sds = sim$source_means_sds, tdf_means_sds = sim$tdf_means_sds,
#'   conc_dep = sim$conc_dep, error_structure = sim$error_structure,
#'   source_col = sim$source_col
#' )
bsimms_get_prior <- function(formula, mixture_data, source_data, tdf_data, isotope_names,
                              source_means_sds = FALSE, tdf_means_sds = TRUE,
                              conc_dep = FALSE,
                              error_structure = c("process_residual", "process_only", "residual_only"),
                              source_col = "Source") {
  spec <- build_bsimms_spec(
    formula = formula, mixture_data = mixture_data,
    source_data = source_data, tdf_data = tdf_data, isotope_names = isotope_names,
    source_means_sds = source_means_sds, tdf_means_sds = tdf_means_sds,
    conc_dep = conc_dep, error_structure = error_structure, source_col = source_col
  )
  default_bsimms_prior(spec)
}
