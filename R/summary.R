#' Print a `bsimms` fit
#'
#' Prints a one-screen overview of a fitted model: formula, source/isotope/
#' mixture counts and names, error structure, source/TDF data mode,
#' whether concentration dependence is enabled, group-level terms (if any),
#' and the fitting backend used. For posterior summaries, use
#' [summary.bsimms_fit()].
#'
#' @param x A `bsimms_fit` object (as returned by [bsimm()]).
#' @param ... Currently unused.
#' @return `x`, invisibly.
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
#' fit
#' }
print.bsimms_fit <- function(x, ...) {
  spec <- x$spec
  cat("Bayesian stable isotope mixing model (bsimms)\n")
  cat(" formula:         ", deparse(spec$formula), "\n", sep = "")
  cat(
    " sources (K):     ",
    spec$K,
    " (",
    paste(spec$source_names, collapse = ", "),
    ")\n",
    sep = ""
  )
  cat(
    " isotopes (J):    ",
    spec$J,
    " (",
    paste(spec$isotope_names, collapse = ", "),
    ")\n",
    sep = ""
  )
  cat(" mixtures (N):    ", spec$N, "\n", sep = "")
  cat(" error structure: ", spec$error_structure, "\n", sep = "")
  cat(" source data:     ", spec$source$mode, "\n", sep = "")
  cat(" tdf data:        ", spec$tdf$mode, "\n", sep = "")
  if (spec$has_conc_dep) {
    cat(" concentration dependence: yes\n")
  }
  if (length(spec$re_terms) > 0) {
    cat(" group-level terms:\n")
    for (re in spec$re_terms) {
      cat(
        "   (",
        paste(re$term_names, collapse = " + "),
        " | ",
        re$group,
        "), ",
        length(re$group_levels),
        " levels\n",
        sep = ""
      )
    }
  }
  cat(" backend:         ", x$backend, "\n", sep = "")
  invisible(x)
}

#' Summarise a `bsimms` fit
#'
#' Reports posterior summaries for the global (population-average) source
#' proportions (`p_global`, the `Dirichlet`-distributed baseline shared by
#' all mixture samples, as in `MixSIAR`), the ILR-scale fixed-effect covariate
#' slopes (deviations from that baseline, if any), group-level standard
#' deviations (if any), and the error term(s).
#'
#' @param object A `bsimms_fit` object (as returned by [bsimm()]).
#' @param robust Logical; if `FALSE` (default) summarise the central
#'   tendency/spread with `mean`/`sd`, if `TRUE` use `median`/`mad` instead.
#' @param probs Quantiles to report. Default 2.5%/97.5%.
#' @param ... Currently unused.
#' @return An object of class `summary.bsimms_fit`.
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
#' summary(fit)
#' }
summary.bsimms_fit <- function(
  object,
  robust = FALSE,
  probs = c(0.025, 0.975),
  ...
) {
  spec <- object$spec
  measures <- c(
    if (robust) list("median", "mad") else list("mean", "sd"),
    list(
      ~ posterior::quantile2(.x, probs = probs),
      "rhat",
      "ess_bulk",
      "ess_tail"
    )
  )
  summarise <- function(draws) {
    as.data.frame(do.call(posterior::summarise_draws, c(list(draws), measures)))
  }

  fixed <- NULL
  if (spec$P > 0) {
    fixed <- summarise(bsimms_draws(object, variable = "beta"))
    fixed$fixed_effect <- rep(spec$fixed_names, spec$D)
    fixed$ilr_dim <- rep(seq_len(spec$D), each = spec$P)
  }

  re_summary <- NULL
  if (length(spec$re_terms) > 0) {
    re_summary <- do.call(
      rbind,
      lapply(spec$re_terms, function(re) {
        d <- summarise(bsimms_draws(
          object,
          variable = paste0("sd_re_", re$label)
        ))
        d$group <- re$group
        d$term <- rep(re$term_names, each = spec$D)
        d$ilr_dim <- rep(seq_len(spec$D), times = length(re$term_names))
        d$variable <- NULL
        d
      })
    )
  }

  err_summary <- NULL
  if (spec$error_structure %in% c("residual_only", "process_residual")) {
    err_var <- if (spec$error_structure == "residual_only") {
      "sigma"
    } else {
      "resid_prop"
    }
    err_summary <- summarise(bsimms_draws(object, variable = err_var))
    err_summary$isotope <- spec$isotope_names
  }

  pop_summary <- summarise(bsimms_draws(object, variable = "p_global"))
  pop_summary$source <- spec$source_names
  pop_summary$variable <- NULL

  out <- list(
    call = object$call,
    spec = spec,
    robust = robust,
    fixed = fixed,
    group = re_summary,
    error = err_summary,
    population_proportions = pop_summary
  )
  class(out) <- "summary.bsimms_fit"
  out
}

#' @export
print.summary.bsimms_fit <- function(x, digits = 3, ...) {
  cat("Bayesian stable isotope mixing model\n")
  cat(" formula: ", deparse(x$spec$formula), "\n\n", sep = "")

  est_cols <- if (isTRUE(x$robust)) c("median", "mad") else c("mean", "sd")
  diag_cols <- c("rhat", "ess_bulk", "ess_tail")
  non_quant_cols <- c(
    "variable",
    "mean",
    "sd",
    "median",
    "mad",
    diag_cols,
    "source",
    "fixed_effect",
    "ilr_dim",
    "isotope",
    "group",
    "term"
  )
  quant_cols <- function(d) setdiff(names(d), non_quant_cols)
  round_table <- function(d, id_cols) {
    for (col in setdiff(names(d), id_cols)) {
      d[[col]] <- round(
        d[[col]],
        if (col %in% c("ess_bulk", "ess_tail")) 0 else digits
      )
    }
    d
  }

  cat("Population-average source proportions:\n")
  pp <- x$population_proportions[, c(
    "source",
    est_cols,
    quant_cols(x$population_proportions),
    diag_cols
  )]
  pp <- round_table(pp, "source")
  print(pp, row.names = FALSE)

  if (!is.null(x$fixed)) {
    cat(
      "\nFixed effects (ILR scale, deviations from the global proportions):\n"
    )
    fe <- x$fixed[, c(
      "fixed_effect",
      "ilr_dim",
      est_cols,
      quant_cols(x$fixed),
      diag_cols
    )]
    names(fe)[1:2] <- c("coefficient", "ilr_dim")
    fe <- round_table(fe, c("coefficient", "ilr_dim"))
    print(fe, row.names = FALSE)
  }

  if (!is.null(x$group)) {
    cat("\nGroup-level standard deviations (ILR scale):\n")
    gg <- x$group[, c(
      "group",
      "term",
      "ilr_dim",
      est_cols,
      quant_cols(x$group),
      diag_cols
    )]
    gg <- round_table(gg, c("group", "term", "ilr_dim"))
    print(gg, row.names = FALSE)
  }

  if (!is.null(x$error)) {
    cat("\nError term(s):\n")
    ee <- x$error[, c("isotope", est_cols, quant_cols(x$error), diag_cols)]
    ee <- round_table(ee, "isotope")
    print(ee, row.names = FALSE)
  }
  invisible(x)
}
