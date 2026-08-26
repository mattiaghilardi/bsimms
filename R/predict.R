#' Posterior source proportions
#'
#' Returns raw posterior draws of source proportions `p`, either for the
#' mixture samples in the fitted data (`newdata = NULL`) or for new mixture
#' covariate combinations. Use [fitted_proportions()] for a summarised
#' (mean/sd or median/mad, with a credible interval) version instead.
#'
#' @param object A `bsimms_fit` object (as returned by [bsimm()]).
#' @param newdata Optional data frame of new mixture covariate values (same
#'   columns as used in `formula`).
#' @param re_formula Which group-level (random-effect) terms to condition on,
#'   whether predicting for the fitted mixture samples (`newdata = NULL`) or
#'   for new ones: `NULL` (default) conditions on every group-level term in
#'   the fitted model (with `newdata`, every term's grouping column(s) must
#'   then be supplied); `NA` or `~0` conditions on none of them
#'   (population-average for every term, regardless of what `newdata`
#'   contains); a reduced formula naming a subset, e.g. `~ (1 | Site)`,
#'   conditions only on the named term(s) (with `newdata`, only those
#'   term(s)' columns are needed). When `newdata` is `NULL`,
#'   `allow_new_levels`/`sample_new_levels` are ignored, since a fitted
#'   mixture sample's group membership is never new.
#' @param allow_new_levels Logical; only relevant with `newdata`. If `FALSE`
#'   (default), a grouping column value not seen when fitting (including
#'   `NA`, which is always treated as a new level) raises an error. If
#'   `TRUE`, a posterior draw is instead sampled for that new level via
#'   `sample_new_levels`; every row sharing the same new, non-`NA` level
#'   value gets the same sampled draws (as for a real level), while each
#'   `NA` row is sampled independently (it asserts no shared identity).
#' @param sample_new_levels How to sample a new level's group-level deviation
#'   when `allow_new_levels = TRUE`: `"uncertainty"` (default) draws, for
#'   each posterior draw, from a randomly chosen *existing* level's draw at
#'   that same iteration -- most appropriate with many existing levels,
#'   where this empirically reflects the observed between-level spread;
#'   `"gaussian"` instead draws a fresh value for the new level from a
#'   normal distribution centred at zero, using that draw's own estimated
#'   group-level SD (and correlation, for multi-column terms) -- more
#'   appropriate with few existing levels, where `"uncertainty"` could only
#'   ever resample one of a handful of specific levels rather than
#'   represent a genuinely new one.
#' @param ndraws Number of posterior draws to use, randomly subset from the
#'   full posterior. `NULL` (default) uses all draws.
#' @param ... Currently unused.
#' @return A numeric `[n_draws, n_obs, K]` array of proportion draws, with
#'   source names attached as the 3rd dimension's `dimnames` (e.g.
#'   `p_arr[, , "Deer"]`).
#' @export
posterior_proportions <- function(object, ...) {
  UseMethod("posterior_proportions")
}

#' @rdname posterior_proportions
#' @export
posterior_proportions.bsimms_fit <- function(object, newdata = NULL, re_formula = NULL,
                                              allow_new_levels = FALSE,
                                              sample_new_levels = c("uncertainty", "gaussian"),
                                              ndraws = NULL, ...) {
  spec <- object$spec
  p_arr <- if (!is.null(newdata)) {
    sample_new_levels <- rlang::arg_match(sample_new_levels)
    predict_p_newdata(object, newdata, re_formula, allow_new_levels, sample_new_levels, ndraws)
  } else if (is.null(re_formula)) {
    dm <- draws_matrix(object, variable = "p")
    dm <- subset_ndraws(dm, ndraws)
    extract_array_draws(dm, "p", spec$N, spec$K)
  } else {
    dm <- draws_matrix(object)
    dm <- subset_ndraws(dm, ndraws)
    candidate_terms <- select_re_terms(spec$re_terms, re_formula)
    compute_p_fitted(spec, dm, candidate_terms)
  }
  dimnames(p_arr) <- list(NULL, NULL, spec$source_names)
  p_arr
}

#' Posterior source proportions (summarised)
#'
#' Returns posterior summaries of source proportions `p`, via
#' [posterior_proportions()], either for the mixture samples in the fitted
#' data (`newdata = NULL`) or for new mixture covariate combinations.
#'
#' @inheritParams posterior_proportions
#' @param summary Logical; return posterior summaries (default) or raw
#'   draws (`FALSE`, equivalent to [posterior_proportions()]).
#' @param robust Logical; if `FALSE` (default) summarise the central
#'   tendency/spread with `mean`/`sd`, if `TRUE` use `median`/`mad` instead.
#' @param probs Quantiles to report when `summary = TRUE`.
#' @return If `summary = TRUE`, a long-format data frame with one row per
#'   (observation, source): columns `row`, `source`, and one column per
#'   summary measure, named as in [posterior::summarise_draws()] (`mean`,
#'   `sd`, `q2.5`, `q97.5`, or with `robust = TRUE`, `median`, `mad`,
#'   `q2.5`, `q97.5`). If `summary = FALSE`, a numeric array
#'   `[n_draws, n_obs, K]`.
#' @export
fitted_proportions <- function(object, newdata = NULL, re_formula = NULL, allow_new_levels = FALSE,
                                sample_new_levels = c("uncertainty", "gaussian"), ndraws = NULL,
                                summary = TRUE, robust = FALSE, probs = c(0.025, 0.975)) {
  spec <- object$spec
  p_arr <- posterior_proportions(
    object,
    newdata = newdata, re_formula = re_formula, allow_new_levels = allow_new_levels,
    sample_new_levels = sample_new_levels, ndraws = ndraws
  )
  if (!summary) return(p_arr)
  summarise_draws_by_row(p_arr, spec$source_names, probs, "source", robust)
}

#' Expected mixture isotope values
#'
#' Returns raw posterior draws of the model's expected mixture isotope
#' value `mu` (the mixing-model prediction *before* observation error is
#' added), either for the fitted mixture samples (`newdata = NULL`) or for
#' new mixture covariate combinations. Use [fitted.bsimms_fit()] for a
#' summarised version instead, [posterior_proportions()] for the underlying
#' source proportions, or [posterior_predict.bsimms_fit()] for posterior
#' predictive draws that include observation error.
#'
#' @inheritParams posterior_proportions
#' @param resp Character vector of isotope name(s) (from `isotope_names`)
#'   to return. `NULL` (default) returns all isotopes.
#' @param ... Currently unused.
#' @return A numeric `[n_draws, n_obs, length(resp)]` array.
#' @exportS3Method rstantools::posterior_epred
posterior_epred.bsimms_fit <- function(object, newdata = NULL, resp = NULL, re_formula = NULL,
                                        allow_new_levels = FALSE,
                                        sample_new_levels = c("uncertainty", "gaussian"),
                                        ndraws = NULL, ...) {
  spec <- object$spec
  resp <- if (is.null(resp)) {
    spec$isotope_names
  } else {
    rlang::arg_match(resp, values = spec$isotope_names, multiple = TRUE)
  }
  j_idx <- match(resp, spec$isotope_names)
  if (!is.null(newdata)) {
    sample_new_levels <- rlang::arg_match(sample_new_levels)
    mu_arr <- predict_mu_newdata(object, newdata, re_formula, allow_new_levels, sample_new_levels, ndraws)
    return(mu_arr[, , j_idx, drop = FALSE])
  }
  if (is.null(re_formula)) {
    dm <- draws_matrix(object, variable = "mu")
    dm <- subset_ndraws(dm, ndraws)
    mu_arr <- extract_array_draws(dm, "mu", spec$N, spec$J)
    return(mu_arr[, , j_idx, drop = FALSE])
  }
  needs <- bsimms_needs_flags(spec)
  dm <- draws_matrix(object)
  dm <- subset_ndraws(dm, ndraws)
  candidate_terms <- select_re_terms(spec$re_terms, re_formula)
  p_arr <- compute_p_fitted(spec, dm, candidate_terms)
  mu_arr <- compute_mu_newdata(spec, dm, p_arr, needs$proc)$mu
  mu_arr[, , j_idx, drop = FALSE]
}

#' Fitted mixture isotope values (summarised)
#'
#' Returns posterior summaries of the model's expected mixture isotope
#' value `mu`, via [posterior_epred.bsimms_fit()].
#'
#' @inheritParams posterior_epred.bsimms_fit
#' @param summary Logical; return posterior summaries (default) or raw
#'   draws (`FALSE`, equivalent to [posterior_epred.bsimms_fit()]).
#' @param robust Logical; if `FALSE` (default) summarise the central
#'   tendency/spread with `mean`/`sd`, if `TRUE` use `median`/`mad` instead.
#' @param probs Quantiles to report when `summary = TRUE`.
#' @return If `summary = TRUE`, a long-format data frame with one row per
#'   (observation, isotope): columns `row`, `isotope`, and one column per
#'   summary measure, named as in [posterior::summarise_draws()] (`mean`,
#'   `sd`, `q2.5`, `q97.5`, or with `robust = TRUE`, `median`, `mad`,
#'   `q2.5`, `q97.5`). If `summary = FALSE`, a numeric array
#'   `[n_draws, n_obs, length(resp)]`.
#' @export
fitted.bsimms_fit <- function(object, newdata = NULL, resp = NULL, re_formula = NULL,
                               allow_new_levels = FALSE,
                               sample_new_levels = c("uncertainty", "gaussian"), ndraws = NULL,
                               summary = TRUE, robust = FALSE, probs = c(0.025, 0.975), ...) {
  spec <- object$spec
  resp <- if (is.null(resp)) {
    spec$isotope_names
  } else {
    rlang::arg_match(resp, values = spec$isotope_names, multiple = TRUE)
  }
  mu_arr <- posterior_epred.bsimms_fit(
    object,
    newdata = newdata, resp = resp, re_formula = re_formula, allow_new_levels = allow_new_levels,
    sample_new_levels = sample_new_levels, ndraws = ndraws
  )
  if (!summary) return(mu_arr)
  summarise_draws_by_row(mu_arr, resp, probs, "isotope", robust)
}

#' Summarise an `[n_draws, n_obs, n_var]` draws array via
#' [posterior::summarise_draws()], one row per `(obs, var)` pair.
#'
#' @param arr Numeric `[n_draws, n_obs, n_var]` array of draws.
#' @param var_names Character vector of variable names (length `n_var`),
#'   e.g. isotope or source names.
#' @param probs Quantiles to report, passed to [posterior::quantile2()].
#' @param var_col Name to give the variable-identity column in the output
#'   (e.g. `"source"` or `"isotope"`).
#' @param robust Logical; if `FALSE` (default) summarise the central
#'   tendency/spread with `mean`/`sd`, if `TRUE` use `median`/`mad` instead.
#' @return A long-format data frame with one row per `(obs, var)` pair:
#'   columns `row`, `<var_col>`, and one column per summary measure, named
#'   as in [posterior::summarise_draws()].
#' @noRd
summarise_draws_by_row <- function(arr, var_names, probs, var_col, robust = FALSE) {
  n_draws <- dim(arr)[1]
  n_obs <- dim(arr)[2]
  measures <- c(
    if (robust) list("median", "mad") else list("mean", "sd"),
    list(~ posterior::quantile2(.x, probs = probs))
  )
  out <- do.call(rbind, lapply(seq_len(n_obs), function(i) {
    d <- matrix(arr[, i, ], nrow = n_draws, dimnames = list(NULL, var_names))
    summ <- as.data.frame(do.call(
      posterior::summarise_draws, c(list(posterior::as_draws_matrix(d)), measures)
    ))
    names(summ)[names(summ) == "variable"] <- var_col
    cbind(data.frame(row = i), summ)
  }))
  rownames(out) <- NULL
  out
}

#' Evaluate a group-level term's grouping value for each row of `newdata`,
#' robustly handling crossed/nested groups (e.g. `"Individual:Site"`, as
#' produced by `(1 | Individual:Site)` or an expanded `(1 | Site/Individual)`)
#' by splitting on `:` and pasting the constituent columns together, rather
#' than evaluating `:` as an R expression. Base R's `:` operator only
#' computes the factor interaction when both operands are already factors
#' (see `parse_bsimms_formula()`), so evaluating it directly against `newdata`
#' silently produces the wrong (or `NA`) values whenever a constituent column
#' is a plain character vector -- the common case for user-supplied `newdata`.
#' Assumes every constituent column is already present (see
#' `check_re_newdata_vars()`), so unlike an earlier version of this
#' function it does not itself handle missing columns.
#'
#' @param group A group-level term's `group` string (`re$group`), e.g.
#'   `"Region"` or `"Individual:Site"`.
#' @param newdata Data frame of new mixture covariate values.
#' @return A character vector (length `nrow(newdata)`), in the same
#'   `"<level1>:<level2>"` format as `re$group_levels`. `NA` in any
#'   constituent column propagates to `NA` in the result (handled by
#'   `resolve_re_b()` as a new level).
#' @noRd
combine_group_var <- function(group, newdata) {
  vars <- strsplit(group, ":")[[1]]
  Reduce(
    function(a, b) paste(a, as.character(newdata[[b]]), sep = ":"),
    vars[-1],
    as.character(newdata[[vars[1]]])
  )
}

#' Resolve `re_formula` into the subset of `re_terms` to condition on.
#'
#' @param re_terms `spec$re_terms` (see `parse_bsimms_formula()`).
#' @param re_formula `NULL` (all terms), `NA`/`~0` (no terms), or a
#'   reduced bars-formula naming a subset (see [posterior_proportions()]).
#' @return The selected subset of `re_terms` (possibly empty).
#' @noRd
select_re_terms <- function(re_terms, re_formula) {
  if (is.null(re_formula)) {
    return(re_terms)
  }
  if (!inherits(re_formula, "formula")) {
    if (length(re_formula) == 1 && is.na(re_formula)) {
      return(list())
    }
    cli::cli_abort(
      "{.arg re_formula} must be `NULL`, `NA`, `~0`, or a formula naming group-level term(s), e.g. {.code ~ (1 | Region)}.",
      call = NULL
    )
  }
  if (identical(deparse(re_formula), "~0")) {
    return(list())
  }
  bars <- reformulas::findbars(reformulas::expandDoubleVerts(re_formula))
  if (is.null(bars)) {
    return(list())
  }
  wanted <- vapply(bars, function(b) trimws(strsplit(rlang::expr_deparse(b), "\\|")[[1]][2]), character(1))
  known <- vapply(re_terms, function(re) re$group, character(1))
  unknown <- setdiff(wanted, known)
  if (length(unknown) > 0) {
    cli::cli_abort(
      "{.arg re_formula} refers to group-level term{?s} not in the fitted model: {.field {unknown}}.",
      call = NULL
    )
  }
  Filter(function(re) re$group %in% wanted, re_terms)
}

#' Validate that every constituent grouping variable of every candidate
#' group-level term (as selected by `re_formula`, see `select_re_terms()`)
#' is present in `newdata`. Since `re_formula = NULL` conditions on every
#' term, this means `newdata` must supply every grouping variable in the
#' model; a reduced `re_formula` only requires the columns its selected
#' term(s) need.
#'
#' @param candidate_terms The subset of `spec$re_terms` selected by
#'   `re_formula` (see `select_re_terms()`).
#' @param newdata Data frame of new mixture covariate values.
#' @return Invisible `NULL`; called for its error side effect.
#' @noRd
check_re_newdata_vars <- function(candidate_terms, newdata) {
  for (re in candidate_terms) {
    vars <- strsplit(re$group, ":")[[1]]
    missing_vars <- setdiff(vars, names(newdata))
    if (length(missing_vars) == 0) next
    if (length(vars) == 1) {
      cli::cli_abort(
        "{.arg newdata} is missing column {.field {re$group}}, needed for its group-level term.",
        call = NULL
      )
    }
    cli::cli_abort(
      c(
        "{.arg newdata} is missing {cli::qty(missing_vars)} column{?s} needed for group-level term {.val {re$group}}:",
        "x" = "{.field {missing_vars}}."
      ),
      call = NULL
    )
  }
}

#' Resolve a group-level term's per-row, per-draw contribution
#' (`[n_draws, size, nrow(newdata)]`), looking up each row's grouping
#' value among the term's fitted levels (`re$group_levels`) or, for a new
#' level (unseen, or `NA` -- always treated as new), erroring
#' (`allow_new_levels = FALSE`) or sampling a value via
#' `sample_new_levels` (`allow_new_levels = TRUE`). Every row sharing the
#' same new, non-`NA` value gets the same sampled draws (matching how a
#' real level's draws are shared across its rows); each `NA` row is
#' sampled independently, since `NA` asserts no shared identity.
#'
#' @param re One element of `spec$re_terms`.
#' @param gvar The term's grouping value per row of `newdata`, as returned
#'   by `combine_group_var()`.
#' @param dm A `posterior::draws_matrix` (see `draws_matrix()`), for the
#'   whole fit (not yet subset to this term's variables).
#' @param n_draws Number of posterior draws (`nrow(dm)`).
#' @param D Number of ILR dimensions (`spec$D`).
#' @param allow_new_levels,sample_new_levels See [posterior_proportions()].
#' @return Numeric `[n_draws, size, nrow(newdata)]` array, `size =
#'   length(re$term_names) * D`.
#' @noRd
resolve_re_b <- function(re, gvar, dm, n_draws, D, allow_new_levels, sample_new_levels) {
  size <- length(re$term_names) * D
  n_group_levels <- length(re$group_levels)
  b_arr <- extract_array_draws(dm, paste0("b_re_", re$label), size, n_group_levels) # n_draws x size x G
  lvl <- match(gvar, re$group_levels)
  is_new <- is.na(lvl)

  b_use <- array(NA_real_, dim = c(n_draws, size, length(gvar)))
  for (i in which(!is_new)) b_use[, , i] <- b_arr[, , lvl[i]]
  if (!any(is_new)) {
    return(b_use)
  }

  if (!allow_new_levels) {
    unseen <- unique(gvar[is_new])
    vars <- strsplit(re$group, ":")[[1]]
    if (length(vars) == 1) {
      cli::cli_abort(
        c(
          "{.arg newdata} contains levels of {.field {re$group}} not seen when fitting the model: {.val {unseen}}.",
          "i" = "Set {.code allow_new_levels = TRUE} to predict for new levels, or use {.arg re_formula} to exclude this term."
        ),
        call = NULL
      )
    }
    cli::cli_abort(
      c(
        "{.arg newdata} contains a combination of {.field {vars}} not seen when fitting the model: {.val {unseen}}.",
        "i" = "Set {.code allow_new_levels = TRUE} to predict for new levels, or use {.arg re_formula} to exclude this term."
      ),
      call = NULL
    )
  }

  # Sampling identity: a repeated, genuinely new (non-NA) value asserts the
  # same underlying level, so shares one sampled draw; `NA` asserts no such
  # identity, so each `NA` row gets its own independent sample.
  sampling_id <- gvar
  na_idx <- which(is.na(gvar))
  sampling_id[na_idx] <- paste0(".newdata_row_", na_idx)
  new_rows <- which(is_new)
  unseen_ids <- unique(sampling_id[new_rows])

  if (sample_new_levels == "gaussian") {
    sd_re_arr <- extract_array_draws(dm, paste0("sd_re_", re$label), size) # n_draws x size
    Lcorr_arr <- if (size >= 2) extract_array_draws(dm, paste0("Lcorr_re_", re$label), size, size) else NULL
  }

  for (u in unseen_ids) {
    rows_u <- new_rows[sampling_id[new_rows] == u]
    new_b_u <- matrix(NA_real_, n_draws, size)
    if (sample_new_levels == "gaussian") {
      for (s in seq_len(n_draws)) {
        z <- stats::rnorm(size)
        new_b_u[s, ] <- if (size >= 2) {
          as.vector(diag(sd_re_arr[s, ], nrow = size, ncol = size) %*% matrix(Lcorr_arr[s, , ], size, size) %*% z)
        } else {
          sd_re_arr[s, ] * z
        }
      }
    } else {
      chosen_g <- sample.int(n_group_levels, n_draws, replace = TRUE)
      for (s in seq_len(n_draws)) new_b_u[s, ] <- b_arr[s, , chosen_g[s]]
    }
    for (i in rows_u) b_use[, , i] <- new_b_u
  }
  b_use
}

#' Compute posterior source-proportion draws for new mixture covariate
#' combinations, from an already-built (and possibly `ndraws`-subset)
#' draws matrix: rebuilds the ILR-scale linear predictor from posterior
#' draws (`ilr_global` always, `beta` if `spec$P > 0`, each candidate
#' group-level term's `b_re_*`, resolved per row via `resolve_re_b()`),
#' then inverse-ILR's it onto the source simplex (via [ilr_inv()]). The
#' core of `predict_p_newdata()`, factored out so `predict_mu_newdata()`/
#' `predict_y_rep_newdata()` can compute `p` from the *same* draws subset
#' they use for source/TDF isotope moments (`ndraws` must only be applied
#' once, since it randomly subsets rows of `dm`).
#'
#' @param spec A `bsimms_spec` (see `build_bsimms_spec()`).
#' @param dm A `posterior::draws_matrix`, already subset to the draws to
#'   use (see `subset_ndraws()`).
#' @param newdata Data frame of new mixture covariate values.
#' @param re_formula,allow_new_levels,sample_new_levels See
#'   [posterior_proportions()].
#' @return Numeric `[nrow(dm), nrow(newdata), K]` array of proportion
#'   draws.
#' @noRd
compute_p_newdata <- function(spec, dm, newdata, re_formula, allow_new_levels, sample_new_levels) {
  candidate_terms <- select_re_terms(spec$re_terms, re_formula)
  check_re_newdata_vars(candidate_terms, newdata)
  n_draws <- nrow(dm)

  # Global (population-average) baseline in ILR space (Egozcue et al. 2003,
  # eq. 25) -- same baseline for every row of `newdata`. Stan already
  # computes and saves this as `ilr_global` (transformed parameters block).
  ilr_global_arr <- extract_array_draws(dm, "ilr_global", spec$D) # n_draws x D

  eta_arr <- array(0, dim = c(n_draws, nrow(newdata), spec$D))
  for (s in seq_len(n_draws)) {
    eta_arr[s, , ] <- matrix(ilr_global_arr[s, ], nrow = nrow(newdata), ncol = spec$D, byrow = TRUE)
  }

  if (spec$P > 0) {
    check_newdata_vars(spec$fixed_formula, newdata)
    X_new <- stats::model.matrix(spec$fixed_formula, data = newdata_model_frame(spec$fixed_formula, spec$fixed_frame, newdata))
    X_new <- align_columns(X_new, spec$fixed_names)
    beta_arr <- extract_array_draws(dm, "beta", spec$P, spec$D) # n_draws x P x D
    for (s in seq_len(n_draws)) eta_arr[s, , ] <- eta_arr[s, , ] + X_new %*% beta_arr[s, , ]
  }

  for (re in candidate_terms) {
    gvar <- combine_group_var(re$group, newdata)
    b_use <- resolve_re_b(re, gvar, dm, n_draws, spec$D, allow_new_levels, sample_new_levels) # n_draws x size x nrow(newdata)
    check_newdata_vars(re$term_formula, newdata)
    Zg_new <- stats::model.matrix(re$term_formula, data = newdata)
    Zg_new <- align_columns(Zg_new, re$term_names)
    for (s in seq_len(n_draws)) {
      for (i in seq_len(nrow(newdata))) {
        for (m in seq_along(re$term_names)) {
          for (d in seq_len(spec$D)) {
            eta_arr[s, i, d] <- eta_arr[s, i, d] + Zg_new[i, m] * b_use[s, (m - 1) * spec$D + d, i]
          }
        }
      }
    }
  }

  p_arr <- array(0, dim = c(n_draws, nrow(newdata), spec$K))
  for (s in seq_len(n_draws)) p_arr[s, , ] <- ilr_inv(matrix(eta_arr[s, , ], ncol = spec$D), V = spec$V)
  p_arr
}

#' Compute posterior source-proportion draws for new mixture covariate
#' combinations: builds (and `ndraws`-subsets) the draws matrix, then
#' delegates to `compute_p_newdata()`. Used by
#' [posterior_proportions.bsimms_fit()] for new-data prediction and by
#' `conditional_effects()` for its prediction grid.
#'
#' @param object A `bsimms_fit` object (as returned by [bsimm()]).
#' @param newdata Data frame of new mixture covariate values.
#' @param re_formula,allow_new_levels,sample_new_levels,ndraws See
#'   [posterior_proportions()].
#' @return Numeric `[n_draws, nrow(newdata), K]` array of proportion draws.
#' @noRd
predict_p_newdata <- function(object, newdata, re_formula, allow_new_levels, sample_new_levels, ndraws = NULL) {
  spec <- object$spec
  dm <- draws_matrix(object)
  dm <- subset_ndraws(dm, ndraws)
  compute_p_newdata(spec, dm, newdata, re_formula, allow_new_levels, sample_new_levels)
}

#' Compute posterior source-proportion draws for the fitted mixture
#' samples, conditioning only on `candidate_terms` (as selected by
#' `re_formula`, see `select_re_terms()`) rather than every group-level
#' term. Unlike `compute_p_newdata()`, no new-level handling is needed: a
#' fitted mixture sample's fixed-effect design row (`spec$X`) and each
#' term's design row/group index (`re$Z`/`re$group_idx`) are exactly the
#' ones used when fitting, so they are always valid.
#'
#' @param spec A `bsimms_spec`.
#' @param dm A `posterior::draws_matrix`, already subset to the draws to
#'   use (see `subset_ndraws()`).
#' @param candidate_terms The subset of `spec$re_terms` selected by
#'   `re_formula` (see `select_re_terms()`).
#' @return Numeric `[nrow(dm), spec$N, K]` array of proportion draws.
#' @noRd
compute_p_fitted <- function(spec, dm, candidate_terms) {
  n_draws <- nrow(dm)
  ilr_global_arr <- extract_array_draws(dm, "ilr_global", spec$D) # n_draws x D

  eta_arr <- array(0, dim = c(n_draws, spec$N, spec$D))
  for (s in seq_len(n_draws)) {
    eta_arr[s, , ] <- matrix(ilr_global_arr[s, ], nrow = spec$N, ncol = spec$D, byrow = TRUE)
  }

  if (spec$P > 0) {
    beta_arr <- extract_array_draws(dm, "beta", spec$P, spec$D) # n_draws x P x D
    for (s in seq_len(n_draws)) eta_arr[s, , ] <- eta_arr[s, , ] + spec$X %*% beta_arr[s, , ]
  }

  for (re in candidate_terms) {
    size <- length(re$term_names) * spec$D
    b_arr <- extract_array_draws(dm, paste0("b_re_", re$label), size, length(re$group_levels)) # n_draws x size x G
    for (s in seq_len(n_draws)) {
      for (i in seq_len(spec$N)) {
        g <- re$group_idx[i]
        for (m in seq_along(re$term_names)) {
          for (d in seq_len(spec$D)) {
            eta_arr[s, i, d] <- eta_arr[s, i, d] + re$Z[i, m] * b_arr[s, (m - 1) * spec$D + d, g]
          }
        }
      }
    }
  }

  p_arr <- array(0, dim = c(n_draws, spec$N, spec$K))
  for (s in seq_len(n_draws)) p_arr[s, , ] <- ilr_inv(matrix(eta_arr[s, , ], ncol = spec$D), V = spec$V)
  p_arr
}

#' Rebuild `newdata`'s factor columns to match the levels seen at fit time
#' (via `xlev`, the same mechanism `stats::predict.lm()` and friends use),
#' before `stats::model.matrix()` sees them. Without this, a factor
#' fixed-effect column in `newdata` is dummy-coded using only whatever
#' levels happen to be present in `newdata` itself -- e.g. a single-row
#' `newdata` produces a single-level factor, which `model.matrix()` cannot
#' contrast at all (`"contrasts can be applied only to factors with 2 or
#' more levels"`), and a genuinely new, unseen level would silently be
#' accepted (and produce a nonsensical design matrix) rather than erroring.
#' Also errors explicitly if a column fitted as a factor is supplied in
#' `newdata` as a non-factor (e.g. numeric), which `xlev` alone only
#' warns about, not errors.
#'
#' @param formula A formula (`spec$fixed_formula`).
#' @param fitted_frame The corresponding model frame at fit time
#'   (`spec$fixed_frame`), used only for its factor columns' levels.
#' @param newdata Data frame of new mixture covariate values.
#' @return `newdata`, with its factor columns releveled to match
#'   `fitted_frame`. Errors if a column fitted as a factor is not a factor
#'   or character in `newdata`, or if it contains a level not seen when
#'   fitting.
#' @noRd
newdata_model_frame <- function(formula, fitted_frame, newdata) {
  is_factor <- vapply(fitted_frame, is.factor, logical(1))
  factor_vars <- intersect(names(fitted_frame)[is_factor], names(newdata))
  bad_type <- factor_vars[!vapply(newdata[factor_vars], function(x) is.factor(x) || is.character(x), logical(1))]
  if (length(bad_type) > 0) {
    cli::cli_abort(
      "{.field {bad_type}} must be a factor or character in {.arg newdata}, as it was when fitting the model.",
      call = NULL
    )
  }
  xlev <- lapply(fitted_frame[is_factor], levels)
  for (v in factor_vars) {
    unseen <- setdiff(as.character(newdata[[v]]), xlev[[v]])
    if (length(unseen) > 0) {
      cli::cli_abort(
        c(
          "{.arg newdata} contains levels of {.field {v}} not seen when fitting the model: {.val {unseen}}.",
          "i" = "Valid levels for {.field {v}} are: {.val {xlev[[v]]}}."
        ),
        call = NULL
      )
    }
  }
  stats::model.frame(formula, data = newdata, xlev = xlev)
}

#' Validate that every variable referenced in a formula is present in
#' `newdata`, before building its design matrix. Without this,
#' `stats::model.matrix()` fails with a raw, formula-variable-specific
#' error (e.g. `"object 'x' not found"`) for the common case of a variable
#' missing entirely -- a different, narrower case from what
#' `align_columns()` catches (a present variable whose `newdata` values
#' produce differently-named dummy columns, e.g. from unseen factor
#' levels).
#'
#' @param formula A formula (`spec$fixed_formula` or a group-level term's
#'   `term_formula`).
#' @param newdata Data frame of new mixture covariate values.
#' @return Invisible `NULL`; called for its error side effect.
#' @noRd
check_newdata_vars <- function(formula, newdata) {
  missing_vars <- setdiff(all.vars(formula), names(newdata))
  if (length(missing_vars) > 0) {
    cli::cli_abort(
      "{.arg newdata} is missing column{?s} needed to build the design matrix: {.field {missing_vars}}.",
      call = NULL
    )
  }
}

#' Reorder (and validate the presence of) a newly-built `model.matrix()`'s
#' columns to match a reference set of column names, e.g. so a `newdata`
#' design matrix lines up with the columns a posterior draws array (`beta`,
#' a group-level term's `Z`) was indexed by at fit time.
#'
#' @param X A `model.matrix()`-built design matrix.
#' @param names_ref Character vector of the expected column names, in the
#'   order the output should follow.
#' @return `X`, subset and reordered to `names_ref`'s columns. Errors if
#'   `X` is missing any column in `names_ref`.
#' @noRd
align_columns <- function(X, names_ref) {
  missing_cols <- setdiff(names_ref, colnames(X))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      "{.arg newdata} is missing column{?s} needed to build the design matrix: {.field {missing_cols}}.",
      call = NULL
    )
  }
  X[, names_ref, drop = FALSE]
}

#' Broadcast a fixed (data, not posterior-varying) matrix to every
#' posterior draw, `[n_draws, nrow(mat), ncol(mat)]`, so downstream code
#' can treat "raw" (posterior parameters) and "summary" (fixed data)
#' source/TDF modes uniformly.
#'
#' @param mat A numeric matrix, e.g. `spec$source$mean`.
#' @param n_draws Number of posterior draws.
#' @return Numeric `[n_draws, nrow(mat), ncol(mat)]` array; `arr[s, , ]`
#'   equals `mat` for every `s`.
#' @noRd
broadcast_fixed_matrix <- function(mat, n_draws) {
  arr <- array(0, dim = c(n_draws, nrow(mat), ncol(mat)))
  for (s in seq_len(n_draws)) arr[s, , ] <- mat
  arr
}

#' Per-draw source isotope mean/SD, `[n_draws, K, J]` each: posterior
#' draws if `spec$source$mode == "raw"`, or the fixed summarised mean/SD
#' broadcast to every draw otherwise (`spec$source$mode == "summary"`;
#' mirrors the `source_mean`/`source_sd` alias `stan_transformed_data_lines()`
#' builds from the equivalent fixed data in that case).
#'
#' @param spec A `bsimms_spec`.
#' @param dm A `posterior::draws_matrix`.
#' @param n_draws Number of posterior draws (`nrow(dm)`).
#' @return A list with `mean` and `sd`, each `[n_draws, K, J]`.
#' @noRd
get_source_moments <- function(spec, dm, n_draws) {
  if (spec$source$mode == "raw") {
    list(
      mean = extract_array_draws(dm, "source_mean", spec$K, spec$J),
      sd = extract_array_draws(dm, "source_sd", spec$K, spec$J)
    )
  } else {
    list(
      mean = broadcast_fixed_matrix(spec$source$mean, n_draws),
      sd = broadcast_fixed_matrix(spec$source$sd, n_draws)
    )
  }
}

#' Per-draw TDF isotope mean/SD, `[n_draws, K, J]` each; see
#' `get_source_moments()` (identical logic, for `spec$tdf`).
#'
#' @inheritParams get_source_moments
#' @return A list with `mean` and `sd`, each `[n_draws, K, J]`.
#' @noRd
get_tdf_moments <- function(spec, dm, n_draws) {
  if (spec$tdf$mode == "raw") {
    list(
      mean = extract_array_draws(dm, "tdf_mean", spec$K, spec$J),
      sd = extract_array_draws(dm, "tdf_sd", spec$K, spec$J)
    )
  } else {
    list(
      mean = broadcast_fixed_matrix(spec$tdf$mean, n_draws),
      sd = broadcast_fixed_matrix(spec$tdf$sd, n_draws)
    )
  }
}

#' Compute the expected mixture isotope value `mu` (and, if `needs_proc`,
#' the propagated source/TDF process variance `proc_var`) for new mixture
#' covariate combinations, from `p_arr` (new-data source-proportion draws,
#' see `compute_p_newdata()`) and posterior draws of source/TDF isotope
#' moments. Mirrors the `mu[i, j]`/`proc_var[i, j]` computation in
#' `stan_transformed_parameters_lines()`.
#'
#' @param spec A `bsimms_spec`.
#' @param dm A `posterior::draws_matrix`, already subset to the draws in
#'   `p_arr`.
#' @param p_arr Numeric `[n_draws, n_new, K]` array (see
#'   `compute_p_newdata()`).
#' @param needs_proc Logical; whether to also compute `proc_var` (see
#'   `bsimms_needs_flags()`).
#' @return A list with `mu` (`[n_draws, n_new, J]`), `proc_var`
#'   (`[n_draws, n_new, J]`, or `NULL` if `!needs_proc`), and `source` /
#'   `tdf` (each `list(mean, sd)`, `[n_draws, K, J]`) for reuse by
#'   `predict_y_rep_newdata()`.
#' @noRd
compute_mu_newdata <- function(spec, dm, p_arr, needs_proc) {
  n_draws <- dim(p_arr)[1]
  n_new <- dim(p_arr)[2]
  source <- get_source_moments(spec, dm, n_draws)
  tdf <- get_tdf_moments(spec, dm, n_draws)

  mu_arr <- array(0, dim = c(n_draws, n_new, spec$J))
  proc_var_arr <- if (needs_proc) array(1e-8, dim = c(n_draws, n_new, spec$J)) else NULL

  for (s in seq_len(n_draws)) {
    pmat <- matrix(p_arr[s, , ], nrow = n_new, ncol = spec$K)
    msum <- matrix(source$mean[s, , ], nrow = spec$K, ncol = spec$J) + matrix(tdf$mean[s, , ], nrow = spec$K, ncol = spec$J)
    if (needs_proc) {
      vsum <- matrix(source$sd[s, , ], nrow = spec$K, ncol = spec$J)^2 + matrix(tdf$sd[s, , ], nrow = spec$K, ncol = spec$J)^2
    }
    if (spec$has_conc_dep) {
      for (j in seq_len(spec$J)) {
        denom <- as.vector(pmat %*% spec$conc[, j])
        pk <- (pmat * matrix(spec$conc[, j], n_new, spec$K, byrow = TRUE)) / denom
        mu_arr[s, , j] <- pk %*% msum[, j]
        if (needs_proc) proc_var_arr[s, , j] <- proc_var_arr[s, , j] + as.vector((pk^2) %*% vsum[, j])
      }
    } else {
      mu_arr[s, , ] <- pmat %*% msum
      if (needs_proc) proc_var_arr[s, , ] <- proc_var_arr[s, , ] + (pmat^2) %*% vsum
    }
  }

  list(mu = mu_arr, proc_var = proc_var_arr, source = source, tdf = tdf)
}

#' Per-draw, per-source isotope covariance `source_cov[k] = L_source_cov[k]
#' %*% t(L_source_cov[k])`, `L_source_cov[k] = diag(source_sd[k, ]) %*%
#' L_source_corr[k]`, `[n_draws, K, J, J]`; only meaningful (and only
#' called) when `needs$mv` (raw source data, 2+ isotopes, source/TDF
#' variance propagated into the mixture). Mirrors the `source_cov`/
#' `L_source_cov` computation in `stan_transformed_parameters_lines()`.
#'
#' @param spec A `bsimms_spec`.
#' @param dm A `posterior::draws_matrix`.
#' @param source_sd_arr Per-draw source isotope SD, `[n_draws, K, J]` (see
#'   `get_source_moments()`).
#' @param n_draws Number of posterior draws (`nrow(dm)`).
#' @return Numeric `[n_draws, K, J, J]` array of per-source covariance
#'   matrices.
#' @noRd
compute_source_cov_newdata <- function(spec, dm, source_sd_arr, n_draws) {
  Lcorr_arr <- extract_array_of_matrices(dm, "L_source_corr", spec$K, spec$J)
  cov_arr <- array(0, dim = c(n_draws, spec$K, spec$J, spec$J))
  for (s in seq_len(n_draws)) {
    for (k in seq_len(spec$K)) {
      L <- diag(source_sd_arr[s, k, ], nrow = spec$J, ncol = spec$J) %*% matrix(Lcorr_arr[s, k, , ], spec$J, spec$J)
      cov_arr[s, k, , ] <- L %*% t(L)
    }
  }
  cov_arr
}

#' Full per-(draw, new row) source/TDF covariance `Omega`, propagated into
#' the mixture: diagonal `proc_var`, off-diagonal cross-isotope terms from
#' `source_cov` (see `compute_source_cov_newdata()`), concentration-weighted
#' if `spec$has_conc_dep`. Mirrors the `Omega` computation in
#' `stan_transformed_parameters_lines()`; only called when `needs$mv`
#' (which implies `spec$J > 1`).
#'
#' @param spec A `bsimms_spec`.
#' @param p_arr Numeric `[n_draws, n_new, K]` array (see
#'   `compute_p_newdata()`).
#' @param proc_var_arr Numeric `[n_draws, n_new, J]` array (see
#'   `compute_mu_newdata()`).
#' @param source_cov_arr Numeric `[n_draws, K, J, J]` array (see
#'   `compute_source_cov_newdata()`).
#' @return Numeric `[n_draws, n_new, J, J]` array of covariance matrices.
#' @noRd
compute_omega_newdata <- function(spec, p_arr, proc_var_arr, source_cov_arr) {
  n_draws <- dim(p_arr)[1]
  n_new <- dim(p_arr)[2]
  omega <- array(0, dim = c(n_draws, n_new, spec$J, spec$J))
  for (s in seq_len(n_draws)) {
    pmat <- matrix(p_arr[s, , ], nrow = n_new, ncol = spec$K)
    for (i in seq_len(n_new)) {
      Om <- diag(proc_var_arr[s, i, ], nrow = spec$J, ncol = spec$J)
      for (j1 in seq_len(spec$J - 1)) {
        for (j2 in (j1 + 1):spec$J) {
          if (spec$has_conc_dep) {
            denom1 <- sum(pmat[i, ] * spec$conc[, j1])
            denom2 <- sum(pmat[i, ] * spec$conc[, j2])
            pk1 <- pmat[i, ] * spec$conc[, j1] / denom1
            pk2 <- pmat[i, ] * spec$conc[, j2] / denom2
            cij <- sum(pk1 * pk2 * source_cov_arr[s, , j1, j2])
          } else {
            cij <- sum(pmat[i, ]^2 * source_cov_arr[s, , j1, j2])
          }
          Om[j1, j2] <- cij
          Om[j2, j1] <- cij
        }
      }
      omega[s, i, , ] <- Om
    }
  }
  omega
}

#' Sample posterior predictive `y_rep` draws for new mixture covariate
#' combinations, from posterior draws of `mu`/`proc_var`/error-structure
#' parameters, matching whichever mixture-likelihood form `needs` selects
#' (see `bsimms_needs_flags()`) -- the same branching
#' `stan_generated_quantities_lines()` uses for the fitted-data `y_rep`.
#' The multivariate cases (`needs$mv`, `needs$resid_cor`) sample via a
#' Cholesky factor of the target covariance (`chol()`, or the fitted
#' model's own `Lcorr_resid` draws) rather than reproducing Stan's exact
#' factor, since only the resulting covariance -- not the specific factor
#' -- matters for sampling.
#'
#' @param spec A `bsimms_spec`.
#' @param dm A `posterior::draws_matrix`, already subset to the draws in
#'   `mu`/`proc_var`.
#' @param p_arr Numeric `[n_draws, n_new, K]` array (see
#'   `compute_p_newdata()`).
#' @param mu_arr,proc_var_arr,source_sd_arr `mu`, `proc_var` (or `NULL`)
#'   and source SD (`[n_draws, K, J]`), as returned by
#'   `compute_mu_newdata()`.
#' @param needs See `bsimms_needs_flags()`.
#' @return Numeric `[n_draws, n_new, J]` array of posterior predictive
#'   draws.
#' @noRd
sample_y_rep_newdata <- function(spec, dm, p_arr, mu_arr, proc_var_arr, source_sd_arr, needs) {
  n_draws <- dim(mu_arr)[1]
  n_new <- dim(mu_arr)[2]
  yrep <- array(0, dim = c(n_draws, n_new, spec$J))

  if (needs$mv) {
    source_cov_arr <- compute_source_cov_newdata(spec, dm, source_sd_arr, n_draws)
    omega <- compute_omega_newdata(spec, p_arr, proc_var_arr, source_cov_arr)
    resid_prop_arr <- if (needs$resid_prop) extract_array_draws(dm, "resid_prop", spec$J) else NULL
    for (s in seq_len(n_draws)) {
      for (i in seq_len(n_new)) {
        Sigma <- matrix(omega[s, i, , ], spec$J, spec$J)
        if (needs$resid_prop) {
          rp_sqrt <- sqrt(resid_prop_arr[s, ])
          Sigma <- Sigma * outer(rp_sqrt, rp_sqrt)
        }
        L <- t(chol(Sigma))
        yrep[s, i, ] <- mu_arr[s, i, ] + as.vector(L %*% stats::rnorm(spec$J))
      }
    }
    return(yrep)
  }

  if (needs$resid_cor) {
    sigma_arr <- extract_array_draws(dm, "sigma", spec$J)
    Lcorr_resid_arr <- extract_array_draws(dm, "Lcorr_resid", spec$J, spec$J)
    for (s in seq_len(n_draws)) {
      L <- diag(sigma_arr[s, ], nrow = spec$J, ncol = spec$J) %*% matrix(Lcorr_resid_arr[s, , ], spec$J, spec$J)
      for (i in seq_len(n_new)) {
        yrep[s, i, ] <- mu_arr[s, i, ] + as.vector(L %*% stats::rnorm(spec$J))
      }
    }
    return(yrep)
  }

  if (needs$proc) {
    resid_prop_arr <- if (needs$resid_prop) extract_array_draws(dm, "resid_prop", spec$J) else NULL
    for (s in seq_len(n_draws)) {
      var_mat <- matrix(proc_var_arr[s, , ], n_new, spec$J)
      if (needs$resid_prop) var_mat <- sweep(var_mat, 2, resid_prop_arr[s, ], `*`)
      mu_vec <- as.vector(matrix(mu_arr[s, , ], n_new, spec$J))
      sd_vec <- as.vector(sqrt(var_mat))
      yrep[s, , ] <- matrix(stats::rnorm(n_new * spec$J, mean = mu_vec, sd = sd_vec), n_new, spec$J)
    }
    return(yrep)
  }

  sigma_arr <- extract_array_draws(dm, "sigma", spec$J)
  for (s in seq_len(n_draws)) {
    mu_vec <- as.vector(matrix(mu_arr[s, , ], n_new, spec$J))
    sd_vec <- rep(sigma_arr[s, ], each = n_new)
    yrep[s, , ] <- matrix(stats::rnorm(n_new * spec$J, mean = mu_vec, sd = sd_vec), n_new, spec$J)
  }
  yrep
}

#' Compute expected mixture isotope value `mu` draws for new mixture
#' covariate combinations: builds (and `ndraws`-subsets) the draws matrix
#' once, computes `p` via `compute_p_newdata()`, then `mu` via
#' `compute_mu_newdata()`. Used by [posterior_epred.bsimms_fit()] for
#' new-data prediction.
#'
#' @param object A `bsimms_fit` object (as returned by [bsimm()]).
#' @param newdata Data frame of new mixture covariate values.
#' @param re_formula,allow_new_levels,sample_new_levels,ndraws See
#'   [posterior_proportions()].
#' @return Numeric `[n_draws, nrow(newdata), J]` array.
#' @noRd
predict_mu_newdata <- function(object, newdata, re_formula, allow_new_levels, sample_new_levels, ndraws = NULL) {
  spec <- object$spec
  needs <- bsimms_needs_flags(spec)
  dm <- draws_matrix(object)
  dm <- subset_ndraws(dm, ndraws)
  p_arr <- compute_p_newdata(spec, dm, newdata, re_formula, allow_new_levels, sample_new_levels)
  compute_mu_newdata(spec, dm, p_arr, needs$proc)$mu
}

#' Sample posterior predictive `y_rep` draws for new mixture covariate
#' combinations: builds (and `ndraws`-subsets) the draws matrix once,
#' computes `p` via `compute_p_newdata()`, `mu`/`proc_var` via
#' `compute_mu_newdata()`, then samples `y_rep` via
#' `sample_y_rep_newdata()`. Used by [posterior_predict.bsimms_fit()] for
#' new-data prediction.
#'
#' @inheritParams predict_mu_newdata
#' @return Numeric `[n_draws, nrow(newdata), J]` array.
#' @noRd
predict_y_rep_newdata <- function(object, newdata, re_formula, allow_new_levels, sample_new_levels, ndraws = NULL) {
  spec <- object$spec
  needs <- bsimms_needs_flags(spec)
  dm <- draws_matrix(object)
  dm <- subset_ndraws(dm, ndraws)
  p_arr <- compute_p_newdata(spec, dm, newdata, re_formula, allow_new_levels, sample_new_levels)
  mu_out <- compute_mu_newdata(spec, dm, p_arr, needs$proc)
  sample_y_rep_newdata(spec, dm, p_arr, mu_out$mu, mu_out$proc_var, mu_out$source$sd, needs)
}

#' Posterior predictive draws of mixture isotope values
#'
#' Returns raw posterior predictive draws of mixture isotope values (i.e.
#' [posterior_epred.bsimms_fit()] plus observation error), either for the
#' fitted mixture samples (`newdata = NULL`) or for new mixture covariate
#' combinations. Use [predict.bsimms_fit()] for a summarised version
#' instead.
#'
#' @inheritParams posterior_epred.bsimms_fit
#' @return A numeric `[n_draws, n_obs, length(resp)]` array.
#' @exportS3Method rstantools::posterior_predict
posterior_predict.bsimms_fit <- function(object, newdata = NULL, resp = NULL, re_formula = NULL,
                                          allow_new_levels = FALSE,
                                          sample_new_levels = c("uncertainty", "gaussian"),
                                          ndraws = NULL, ...) {
  spec <- object$spec
  resp <- if (is.null(resp)) {
    spec$isotope_names
  } else {
    rlang::arg_match(resp, values = spec$isotope_names, multiple = TRUE)
  }
  j_idx <- match(resp, spec$isotope_names)
  if (!is.null(newdata)) {
    sample_new_levels <- rlang::arg_match(sample_new_levels)
    yrep_arr <- predict_y_rep_newdata(object, newdata, re_formula, allow_new_levels, sample_new_levels, ndraws)
    return(yrep_arr[, , j_idx, drop = FALSE])
  }
  if (is.null(re_formula)) {
    dm <- draws_matrix(object, variable = "y_rep")
    dm <- subset_ndraws(dm, ndraws)
    yrep_arr <- extract_array_draws(dm, "y_rep", spec$N, spec$J)
    return(yrep_arr[, , j_idx, drop = FALSE])
  }
  needs <- bsimms_needs_flags(spec)
  dm <- draws_matrix(object)
  dm <- subset_ndraws(dm, ndraws)
  candidate_terms <- select_re_terms(spec$re_terms, re_formula)
  p_arr <- compute_p_fitted(spec, dm, candidate_terms)
  mu_out <- compute_mu_newdata(spec, dm, p_arr, needs$proc)
  yrep_arr <- sample_y_rep_newdata(spec, dm, p_arr, mu_out$mu, mu_out$proc_var, mu_out$source$sd, needs)
  yrep_arr[, , j_idx, drop = FALSE]
}

#' Posterior predictive draws of mixture isotope values (summarised)
#'
#' Returns posterior summaries of the model's posterior predictive draws
#' (i.e. [fitted.bsimms_fit()] plus observation error), via
#' [posterior_predict.bsimms_fit()].
#'
#' @inheritParams posterior_predict.bsimms_fit
#' @param summary Logical; return posterior summaries (default) or raw
#'   draws (`FALSE`, equivalent to [posterior_predict.bsimms_fit()]).
#' @param robust Logical; if `FALSE` (default) summarise the central
#'   tendency/spread with `mean`/`sd`, if `TRUE` use `median`/`mad` instead.
#' @param probs Quantiles to report when `summary = TRUE`.
#' @return If `summary = TRUE`, a long-format data frame (`row`, `isotope`,
#'   and one column per summary measure, named as in
#'   [posterior::summarise_draws()]); if `FALSE`, an array
#'   `[n_draws, n_obs, length(resp)]`.
#' @export
predict.bsimms_fit <- function(object, newdata = NULL, resp = NULL, re_formula = NULL,
                                allow_new_levels = FALSE,
                                sample_new_levels = c("uncertainty", "gaussian"), ndraws = NULL,
                                summary = TRUE, robust = FALSE, probs = c(0.025, 0.975), ...) {
  spec <- object$spec
  resp <- if (is.null(resp)) {
    spec$isotope_names
  } else {
    rlang::arg_match(resp, values = spec$isotope_names, multiple = TRUE)
  }
  yrep_arr <- posterior_predict.bsimms_fit(
    object,
    newdata = newdata, resp = resp, re_formula = re_formula, allow_new_levels = allow_new_levels,
    sample_new_levels = sample_new_levels, ndraws = ndraws
  )
  if (!summary) return(yrep_arr)
  summarise_draws_by_row(yrep_arr, resp, probs, "isotope", robust)
}
