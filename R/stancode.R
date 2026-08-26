#' Generate the Stan code for a `bsimms` model
#'
#' Builds a complete, human-readable Stan program implementing the stable
#' isotope mixing model described by `formula` and the supplied data:
#' nothing is hidden inside compiled internals, the generated `.stan` text
#' is the model, and it can be inspected, hand-edited, and compiled with
#' `cmdstanr` or `rstan` directly.
#'
#' @param formula A one-sided `lme4`-style formula for the source
#'   proportions, e.g. `~ Sex + Season + (1 + Season | Region)`. All
#'   variables referenced must be columns of `mixture_data`. The
#'   population-level baseline is not part of `formula`: it is `p_global`,
#'   a `simplex[K]` parameter with a `Dirichlet(alpha)` prior (as in
#'   `MixSIAR`), and `formula`'s terms are covariate deviations from it,
#'   added in ILR space (see [bsimms_prior()]'s `"p_global"` class to set
#'   `alpha`). A formula with no covariates (e.g. `~ 1`) is a valid
#'   Dirichlet-baseline-only model.
#' @param mixture_data Data frame of mixture observations (e.g. consumers,
#'   in a diet-mixing study): one row per individual (or per sample), with
#'   a column for each entry of `isotope_names`, plus any
#'   covariate/grouping columns used in `formula`.
#' @param source_data Source isotope data. If `source_means_sds = FALSE`
#'   (default), long format with one row per raw sample and columns
#'   `Source`, `isotope_names`. If `source_means_sds = TRUE`, one row per
#'   source with columns `Source`, `<isotope>_mean`, `<isotope>_sd` for
#'   every isotope.
#' @param tdf_data Trophic discrimination factor (diet-tissue discrimination)
#'   data, same layout convention as `source_data`, controlled by
#'   `tdf_means_sds` (default `TRUE`, since TDFs are usually taken from the
#'   literature as means/SDs).
#' @param isotope_names Character vector of isotope column names shared by
#'   `mixture_data`, `source_data` and `tdf_data`, e.g. `c("d13C", "d15N")`.
#' @param source_means_sds Logical; is `source_data` supplied as
#'   means/SDs (`TRUE`) or raw replicate samples (`FALSE`, default)? When
#'   raw, source means and SDs are estimated as part of the model, and
#'   their uncertainty is fully propagated into the posterior.
#' @param tdf_means_sds Logical; is `tdf_data` supplied as means/SDs
#'   (`TRUE`, default) or raw replicate samples (`FALSE`)? When raw, TDF
#'   means and SDs are estimated as part of the model and their uncertainty
#'   is fully propagated into the posterior. Unlike raw source data, though,
#'   no cross-isotope correlation is estimated for raw TDF data (there is
#'   no `"tdf_cor"` class in [bsimms_prior()]): the mixture likelihood only
#'   takes the multivariate (`multi_normal_cholesky`) form when *source*
#'   data is raw and there are 2+ isotopes. This asymmetry is intentional
#'   (raw TDF data with enough replicates per source to usefully estimate a
#'   correlation is a rare combination in practice), not an oversight, but
#'   could be added symmetrically to `source_cor` if needed.
#' @param conc_dep Logical; enable elemental concentration dependence
#'   (default `FALSE`)? If `TRUE`, `source_data` must have a `<isotope>_conc`
#'   column for every isotope in `isotope_names`, giving that source's
#'   proportional elemental concentration (in `(0, 1]`) for that isotope's
#'   element -- one value per raw sample (averaged per source) if
#'   `source_means_sds = FALSE`, or one value per source if
#'   `source_means_sds = TRUE`. Since these are proportions of the same
#'   source's total mass, they cannot sum to more than 1 across isotopes.
#' @param error_structure One of `"process_residual"` (default: source/TDF
#'   variance propagated into the mixture and scaled by an estimated
#'   multiplicative residual-error factor, i.e. the "Residual * Process"
#'   error structure of Stock & Semmens 2016), `"process_only"` (only
#'   source/TDF variance propagated, no separate residual term, as in the
#'   original MixSIR model, Moore & Semmens 2008), or `"residual_only"`
#'   (source/TDF means treated as fixed and known; all unexplained
#'   variance goes into a single residual error term, as in the original
#'   SIAR model, Parnell et al. 2010).
#' @param prior Optional `bsimms_prior` object (see [bsimms_prior()]) with
#'   one or more rows overriding the default priors. Unspecified
#'   parameters keep their (weakly informative, partly data-scaled)
#'   defaults; see [bsimms_get_prior()].
#' @param source_col Name of the source-identifier column shared by
#'   `source_data` and `tdf_data`. Default `"Source"`.
#' @return A single character string of Stan code (class `bsimms_stancode`;
#'   `print()`s as plain text).
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
#'   97(10), 2562-2569. \doi{10.1002/ecy.1517}
#' @export
#' @examples
#' mixture_data <- data.frame(
#'   d13C = c(-20, -21, -19, -22, -20.5, -21.5),
#'   d15N = c(10, 11, 9, 12, 10.5, 11.5),
#'   Region = factor(rep(c("A", "B"), each = 3))
#' )
#' source_data <- data.frame(
#'   Source = rep(c("Beaver", "Deer"), each = 3),
#'   d13C = c(-25, -24, -26, -18, -17, -19),
#'   d15N = c(5, 6, 4, 8, 9, 7)
#' )
#' tdf_data <- data.frame(
#'   Source = c("Beaver", "Deer"),
#'   d13C_mean = c(1, 1.2), d13C_sd = c(0.2, 0.3),
#'   d15N_mean = c(3, 3.1), d15N_sd = c(0.4, 0.5)
#' )
#' code <- make_stancode(
#'   ~ 1 + (1 | Region), mixture_data = mixture_data, source_data = source_data,
#'   tdf_data = tdf_data, isotope_names = c("d13C", "d15N")
#' )
#' cat(code)
make_stancode <- function(formula, mixture_data, source_data, tdf_data, isotope_names,
                           source_means_sds = FALSE, tdf_means_sds = TRUE,
                           conc_dep = FALSE,
                           error_structure = c("process_residual", "process_only", "residual_only"),
                           prior = NULL,
                           source_col = "Source") {
  error_structure <- rlang::arg_match(error_structure)
  spec <- build_bsimms_spec(
    formula = formula, mixture_data = mixture_data,
    source_data = source_data, tdf_data = tdf_data, isotope_names = isotope_names,
    source_means_sds = source_means_sds, tdf_means_sds = tdf_means_sds,
    conc_dep = conc_dep, error_structure = error_structure, source_col = source_col
  )
  prior_df <- merge_bsimms_prior(default_bsimms_prior(spec), prior)
  code <- bsimms_stancode_from_spec(spec, prior_df)
  class(code) <- c("bsimms_stancode", "character")
  code
}

#' Print generated Stan code
#'
#' Prints a `bsimms_stancode` object as plain Stan program text (rather
#' than as a quoted, escaped R character string, which is how it would
#' print without this method).
#'
#' @param x A `bsimms_stancode` object (as returned by [make_stancode()]).
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.bsimms_stancode <- function(x, ...) {
  cat(x)
  invisible(x)
}

#' Derive which optional model components a `bsimms_spec` needs: `sigma`
#' (a single residual-error SD per isotope; `error_structure ==
#' "residual_only"`), `resid_prop` (MixSIAR's `resid.prop` multiplicative
#' factor; `error_structure == "process_residual"`), `proc` (source/TDF
#' variance propagated into a per-mixture-sample process variance;
#' `error_structure %in% c("process_only", "process_residual")`),
#' `source_cor` (raw source data with 2+ isotopes, so per-source isotope
#' correlations are estimated), `mv` (`source_cor && proc`: the mixture
#' needs a full multivariate process + residual covariance rather than
#' per-isotope normals), and `resid_cor` (`sigma` with 2+ isotopes: a
#' shared residual-error correlation is estimated). Shared by
#' `bsimms_stancode_from_spec()` (which of the 7 Stan program blocks'
#' optional pieces to generate) and `predict.R`'s new-data `mu`/`y_rep`
#' prediction (which R computation mirrors which Stan block).
#'
#' @param spec A `bsimms_spec` (see `build_bsimms_spec()`).
#' @return A named list of logicals: `sigma`, `resid_prop`, `proc`,
#'   `source_cor`, `mv`, `resid_cor`.
#' @noRd
bsimms_needs_flags <- function(spec) {
  needs_sigma      <- spec$error_structure == "residual_only"
  needs_resid_prop <- spec$error_structure == "process_residual"
  needs_proc       <- spec$error_structure %in% c("process_only", "process_residual")
  needs_source_cor <- spec$source$mode == "raw" && spec$J > 1
  needs_mv         <- needs_source_cor && needs_proc
  needs_resid_cor  <- needs_sigma && spec$J > 1
  list(
    sigma = needs_sigma, resid_prop = needs_resid_prop, proc = needs_proc,
    source_cor = needs_source_cor, mv = needs_mv, resid_cor = needs_resid_cor
  )
}

#' Assemble the complete Stan program text for a model: computes the
#' `needs_*` flags (see `bsimms_needs_flags()`) that gate which optional
#' pieces of the model are generated, then builds and concatenates each of
#' the 7 Stan program blocks (`functions`, `data`, `transformed data`,
#' `parameters`, `transformed parameters`, `model`, `generated quantities`)
#' via the `stan_*_lines()` helpers below, indented one level with
#' `indent()`.
#'
#' @param spec A `bsimms_spec` (see `build_bsimms_spec()`).
#' @param prior_df A merged `bsimms_prior` data frame (as returned by
#'   `merge_bsimms_prior()`), used by `stan_model_lines()` to look up each
#'   parameter's prior via `select_prior()`.
#' @return A single character string: the complete Stan program.
#' @noRd
bsimms_stancode_from_spec <- function(spec, prior_df) {
  needs <- bsimms_needs_flags(spec)
  needs_sigma      <- needs$sigma
  needs_resid_prop <- needs$resid_prop
  needs_proc       <- needs$proc
  needs_source_cor <- needs$source_cor
  needs_mv         <- needs$mv
  needs_resid_cor  <- needs$resid_cor

  blocks <- c(
    stan_header(spec),
    "functions {", indent(stan_functions_lines(spec)), "}",
    "data {",   indent(stan_data_lines(spec)),   "}",
    "transformed data {", indent(stan_transformed_data_lines(spec)), "}",
    "parameters {", indent(stan_parameters_lines(spec, needs_sigma, needs_resid_prop, needs_source_cor, needs_resid_cor)), "}",
    "transformed parameters {", indent(stan_transformed_parameters_lines(spec, needs_proc, needs_resid_prop, needs_source_cor, needs_mv, needs_resid_cor)), "}",
    "model {", indent(stan_model_lines(spec, prior_df, needs_sigma, needs_resid_prop, needs_proc, needs_source_cor, needs_mv, needs_resid_cor)), "}",
    "generated quantities {", indent(stan_generated_quantities_lines(spec, needs_sigma, needs_resid_prop, needs_proc, needs_mv, needs_resid_cor)), "}"
  )
  paste(blocks, collapse = "\n")
}

#' Prefix each element of a character vector with `n` spaces, for nesting
#' Stan code inside a block (or a block inside another, e.g. a `for` loop
#' body).
#'
#' @param lines Character vector of Stan code lines.
#' @param n Number of spaces to prefix (default 2).
#' @return Character vector, same length as `lines`.
#' @noRd
indent <- function(lines, n = 2) {
  if (length(lines) == 0) return(character(0))
  pad <- strrep(" ", n)
  paste0(pad, lines)
}

#' Build the leading comment-only lines of the Stan program: a
#' human-readable summary of the model (sources, isotopes, fixed/
#' group-level terms, error structure, source/TDF data mode).
#'
#' @param spec A `bsimms_spec` (see `build_bsimms_spec()`).
#' @return Character vector of Stan comment lines (`// ...`).
#' @noRd
stan_header <- function(spec) {
  c(
    paste0("// Stan program generated by bsimms ", utils::packageVersion("bsimms")),
    paste0("// sources: ", paste(spec$source_names, collapse = ", ")),
    paste0("// isotopes: ", paste(spec$isotope_names, collapse = ", ")),
    "// global proportions: p_global ~ Dirichlet(alpha), covariate effects added in ILR space",
    paste0("// fixed effects (covariate slopes, no intercept): ", if (length(spec$fixed_names) > 0) paste(spec$fixed_names, collapse = ", ") else "(none)"),
    if (length(spec$re_terms) > 0) {
      paste0("// group-level terms: ", paste(vapply(spec$re_terms, function(r) {
        sprintf("(%s | %s)", paste(r$term_names, collapse = " + "), r$group)
      }, character(1)), collapse = ", "))
    } else {
      "// group-level terms: (none)"
    },
    paste0("// error structure: ", spec$error_structure),
    paste0("// source data: ", spec$source$mode, if (spec$has_conc_dep) ", concentration-dependent" else ""),
    paste0("// tdf data: ", spec$tdf$mode),
    ""
  )
}

## ---- functions block ----------------------------------------------------

#' Build the `functions` block: the two Stan helper functions used
#' throughout the generated program, `gmean()` (geometric mean of a
#' composition's leading `k` parts, for the forward ILR transform of
#' `p_global`) and `inverse_ilr()` (maps ILR coordinates back onto the
#' source simplex, Egozcue et al. 2003 eq. 24). Currently the same
#' regardless of `spec` (`spec` is unused, kept only so every `stan_*_lines`
#' helper shares a consistent `(spec, ...)` signature).
#'
#' @param spec A `bsimms_spec` (see `build_bsimms_spec()`); currently unused.
#' @return Character vector of Stan code lines.
#' @noRd
stan_functions_lines <- function(spec) {
  c(
    "/* Geometric mean of the first k parts of a composition",
    " * Used for the forward ILR transform of p_global (Egozcue et al. 2003,",
    " * eq. 25)",
    " * Args:",
    " *   p: a composition (simplex-scale vector)",
    " *   k: number of leading parts of p to average",
    " * Returns:",
    " *   the geometric mean of p[1:k]",
    " */",
    "real gmean(vector p, int k) {",
    indent("return exp(mean(log(p[1:k])));"),
    "}",
    "",
    "/* Inverse isometric log-ratio transform",
    " * Maps ILR coordinates back onto the source simplex following Egozcue",
    " * et al. (2003, eq. 24)",
    " * Args:",
    " *   z: ILR-scale coordinates (length K - 1)",
    " *   E: simplex-domain ILR basis, columns e_1, ..., e_{K-1} (K x (K - 1))",
    " * Returns:",
    " *   a simplex of length K",
    " */",
    "vector inverse_ilr(vector z, matrix E) {",
    indent("int K = rows(E);"),
    indent("int D = cols(E);"),
    indent("matrix[K, D] cross;"),
    indent("vector[K] p;"),
    indent("for (d in 1:D) {"),
    indent(indent("vector[K] powered = pow(E[, d], z[d]);")),
    indent(indent("cross[, d] = powered / sum(powered);")),
    indent("}"),
    indent("p = cross[, 1];"),
    indent("if (D > 1) {"),
    indent(indent("for (d in 2:D) {")),
    indent(indent(indent("p = p .* cross[, d];"))),
    indent(indent("}")),
    indent("}"),
    indent("return p / sum(p);"),
    "}"
  )
}

## ---- data block -------------------------------------------------------

#' Build the `data` block: mixture dimensions and data
#' (`N`/`J`/`K`/`D`/`V`/`y`/`alpha_dirichlet`), the fixed-effect design
#' (`P`/`X`, only if `spec$P > 0`), one group of declarations per
#' group-level term in `spec$re_terms`, source and TDF data (raw replicate
#' samples or means/SDs, depending on `spec$source`/`spec$tdf`'s `mode`),
#' and, if `spec$has_conc_dep`, the concentration matrix.
#'
#' @param spec A `bsimms_spec` (see `build_bsimms_spec()`).
#' @return Character vector of Stan code lines.
#' @noRd
stan_data_lines <- function(spec) {
  lines <- c(
    "int<lower=1> N;  // number of mixture observations",
    "int<lower=1> J;  // number of isotopes/tracers",
    "int<lower=2> K;  // number of sources",
    "int<lower=1> D;  // K - 1, ILR dimension",
    "matrix[K, D] V;  // ILR basis matrix",
    "matrix[N, J] y;  // mixture isotope data",
    "vector<lower=0>[K] alpha_dirichlet;  // Dirichlet concentration for p_global"
  )
  if (spec$P > 0) {
    lines <- c(lines,
      "int<lower=1> P;  // number of fixed-effect coefficients (no intercept)",
      "matrix[N, P] X;  // fixed-effect design matrix (no intercept column)"
    )
  }

  for (re in spec$re_terms) {
    Mv <- length(re$term_names)
    Gv <- length(re$group_levels)
    lines <- c(lines,
      sprintf("int<lower=1> N_re_%s;  // levels of grouping factor '%s'", re$label, re$group),
      sprintf("int<lower=1> M_re_%s;  // number of group-level terms for '%s'", re$label, re$group),
      sprintf("matrix[N, M_re_%s] Z_re_%s;  // group-level design matrix for '%s'", re$label, re$label, re$group),
      sprintf("array[N] int<lower=1, upper=N_re_%s> grp_re_%s;  // level index per observation", re$label, re$label)
    )
  }

  if (spec$source$mode == "raw") {
    lines <- c(lines,
      "int<lower=1> N_source_raw;  // number of raw source samples",
      "array[N_source_raw] int<lower=1, upper=K> source_idx;  // source index (1..K) per raw sample",
      "matrix[N_source_raw, J] source_raw;  // raw source isotope measurements"
    )
  } else {
    lines <- c(lines,
      "matrix[K, J] source_mean_data;  // source isotope means, one row per source",
      "matrix[K, J] source_sd_data;  // source isotope SDs, one row per source"
    )
  }

  if (spec$tdf$mode == "raw") {
    lines <- c(lines,
      "int<lower=1> N_tdf_raw;  // number of raw TDF samples",
      "array[N_tdf_raw] int<lower=1, upper=K> tdf_idx;  // source index (1..K) per raw TDF sample",
      "matrix[N_tdf_raw, J] tdf_raw;  // raw TDF (trophic discrimination factor) measurements"
    )
  } else {
    lines <- c(lines,
      "matrix[K, J] tdf_mean_data;  // TDF isotope means, one row per source",
      "matrix[K, J] tdf_sd_data;  // TDF isotope SDs, one row per source"
    )
  }

  if (spec$has_conc_dep) {
    lines <- c(lines, "matrix<lower=0, upper=1>[K, J] conc;  // elemental concentration by source x isotope")
  }

  lines
}

## ---- transformed data --------------------------------------------------

#' Build the `transformed data` block: derives `E`, the simplex-domain ILR
#' basis (`E[, d] = softmax(V[, d])`, Egozcue et al. 2003 eq. 18), from the
#' `V` basis matrix supplied as data; and, for source/TDF data supplied as
#' means/SDs (`spec$source`/`spec$tdf`'s `mode == "summary"`), aliases the
#' `_data`-suffixed input matrices to the plain `source_mean`/`source_sd`/
#' `tdf_mean`/`tdf_sd` names used everywhere else in the program (mirroring
#' the parameters that exist in that role when the corresponding data are
#' raw instead).
#'
#' @param spec A `bsimms_spec` (see `build_bsimms_spec()`).
#' @return Character vector of Stan code lines.
#' @noRd
stan_transformed_data_lines <- function(spec) {
  lines <- c(
    "matrix[K, D] E;  // simplex-domain ILR basis (Egozcue et al. 2003, eq. 18)",
    "for (d in 1:D) {",
    indent("E[, d] = softmax(V[, d]);"),
    "}"
  )
  if (spec$source$mode == "summary") {
    lines <- c(lines, "",
      "matrix[K, J] source_mean = source_mean_data;  // alias: source means as supplied",
      "matrix[K, J] source_sd = source_sd_data;  // alias: source SDs as supplied"
    )
  }
  if (spec$tdf$mode == "summary") {
    lines <- c(lines, "",
      "matrix[K, J] tdf_mean = tdf_mean_data;  // alias: TDF means as supplied",
      "matrix[K, J] tdf_sd = tdf_sd_data;  // alias: TDF SDs as supplied"
    )
  }
  lines
}

## ---- parameters ---------------------------------------------------------

#' Build the `parameters` block: `p_global` (the Dirichlet-distributed
#' population baseline) always; `beta` (fixed-effect slopes) if
#' `spec$P > 0`; each group-level term's non-centered `z_re_*`/`sd_re_*`
#' (and `Lcorr_re_*` if the term has 2+ columns) for every entry in
#' `spec$re_terms`; the error-structure-specific parameters gated by
#' `needs_sigma`/`needs_resid_cor`/`needs_resid_prop`; and, for raw source/
#' TDF data (`spec$source`/`spec$tdf`'s `mode == "raw"`), the estimated
#' `source_mean`/`source_sd`/`tdf_mean`/`tdf_sd` (plus per-source isotope
#' correlation Cholesky factors if `needs_source_cor`).
#'
#' @param spec A `bsimms_spec` (see `build_bsimms_spec()`).
#' @param needs_sigma,needs_resid_prop,needs_source_cor,needs_resid_cor
#'   See `bsimms_stancode_from_spec()`.
#' @return Character vector of Stan code lines.
#' @noRd
stan_parameters_lines <- function(spec, needs_sigma, needs_resid_prop, needs_source_cor, needs_resid_cor) {
  lines <- c("simplex[K] p_global;  // global (population-average) source proportions")
  if (spec$P > 0) {
    lines <- c(lines, "matrix[P, D] beta;  // fixed-effect covariate slopes (no intercept), one column per ILR dimension")
  }

  for (re in spec$re_terms) {
    size <- length(re$term_names) * spec$D
    lines <- c(lines, sprintf("matrix[%d, N_re_%s] z_re_%s;  // std-normal, non-centered", size, re$label, re$label))
    lines <- c(lines, sprintf("vector<lower=0>[%d] sd_re_%s;  // group-level SD, one per (term x ILR dim)", size, re$label))
    if (size > 1) {
      lines <- c(lines, sprintf("cholesky_factor_corr[%d] Lcorr_re_%s;  // group-level correlation (Cholesky factor)", size, re$label))
    }
  }

  if (needs_sigma) {
    lines <- c(lines, "vector<lower=0>[J] sigma;  // residual / observation error")
  }

  if (needs_resid_cor) {
    lines <- c(lines, "cholesky_factor_corr[J] Lcorr_resid;  // residual-error correlation (Cholesky factor)")
  }

  if (needs_resid_prop) {
    lines <- c(lines, "vector<lower=0, upper=20>[J] resid_prop;  // MixSIAR residual-error factor, scales process variance")
  }

  if (spec$source$mode == "raw") {
    lines <- c(lines,
      "matrix[K, J] source_mean;  // estimated source isotope means (raw source data)",
      "matrix<lower=0>[K, J] source_sd;  // estimated source isotope SDs (raw source data)"
    )
    if (needs_source_cor) {
      lines <- c(lines, "array[K] cholesky_factor_corr[J] L_source_corr;  // per-source isotope correlation (Cholesky factor)")
    }
  }
  if (spec$tdf$mode == "raw") {
    lines <- c(lines,
      "matrix[K, J] tdf_mean;  // estimated TDF isotope means (raw TDF data)",
      "matrix<lower=0>[K, J] tdf_sd;  // estimated TDF isotope SDs (raw TDF data)"
    )
  }

  lines
}

## ---- transformed parameters ---------------------------------------------

#' Build the `transformed parameters` block: `ilr_global` (forward ILR of
#' `p_global`, Egozcue et al. 2003 eq. 25) and the per-mixture-sample ILR
#' linear predictor `eta` (baseline, plus `X %*% beta` if `spec$P > 0`,
#' plus each group-level term's contribution from `spec$re_terms`); `p`
#' (`eta` inverse-ILR'd onto the source simplex, one row per mixture
#' sample); `mu` (the expected mixture isotope value, mixing `p` with
#' source/TDF means, optionally concentration-weighted if
#' `spec$has_conc_dep`); and, gated by
#' `needs_proc`/`needs_source_cor`/`needs_mv`/`needs_resid_cor`, the
#' propagated process variance and/or covariance Cholesky factors feeding
#' the mixture likelihood built by `stan_model_lines()`.
#'
#' @param spec A `bsimms_spec` (see `build_bsimms_spec()`).
#' @param needs_proc,needs_resid_prop,needs_source_cor,needs_mv,needs_resid_cor
#'   See `bsimms_stancode_from_spec()`.
#' @return Character vector of Stan code lines.
#' @noRd
stan_transformed_parameters_lines <- function(spec, needs_proc, needs_resid_prop, needs_source_cor, needs_mv, needs_resid_cor) {
  lines <- c(
    "vector[D] ilr_global;  // forward ILR of p_global (Egozcue et al. 2003, eq. 25)",
    "for (d in 1:D) {",
    indent("ilr_global[d] = sqrt(d / (d + 1.0)) * log(gmean(p_global, d) / p_global[d + 1]);"),
    "}",
    "",
    if (spec$P > 0) {
      "matrix[N, D] eta = rep_matrix(to_row_vector(ilr_global), N) + X * beta;  // per-mixture-sample ILR predictor: global baseline + covariate effects"
    } else {
      "matrix[N, D] eta = rep_matrix(to_row_vector(ilr_global), N);  // per-mixture-sample ILR predictor: global baseline only (no covariates)"
    }
  )

  for (re in spec$re_terms) {
    Mv <- length(re$term_names)
    size <- Mv * spec$D
    lines <- c(lines, "", sprintf("// group-level term: (%s | %s)", paste(re$term_names, collapse = " + "), re$group))
    if (size > 1) {
      lines <- c(lines, sprintf("matrix[%d, N_re_%s] b_re_%s = diag_pre_multiply(sd_re_%s, Lcorr_re_%s) * z_re_%s;  // scaled, correlated group-level effects",
                                 size, re$label, re$label, re$label, re$label, re$label))
    } else {
      lines <- c(lines, sprintf("matrix[1, N_re_%s] b_re_%s = rep_matrix(sd_re_%s[1], 1, N_re_%s) .* z_re_%s;  // scaled group-level effects",
                                 re$label, re$label, re$label, re$label, re$label))
    }
    lines <- c(lines,
      sprintf("for (re_i in 1:N) {  // add this term's group-level effect to each mixture sample's eta"),
      indent(sprintf("for (re_m in 1:M_re_%s) {", re$label)),
      indent(indent(sprintf("for (re_d in 1:D) {"))),
      indent(indent(indent(sprintf(
        "eta[re_i, re_d] += Z_re_%s[re_i, re_m] * b_re_%s[(re_m - 1) * D + re_d, grp_re_%s[re_i]];",
        re$label, re$label, re$label
      )))),
      indent(indent("}")),
      indent("}"),
      "}"
    )
  }

  lines <- c(lines, "",
    "array[N] simplex[K] p;  // source proportions, one simplex per mixture sample",
    "for (i in 1:N) {  // inverse-ILR each mixture sample's eta onto the source simplex",
    indent("p[i] = inverse_ilr(to_vector(eta[i]), E);  // Egozcue et al. 2003, eq. 24"),
    "}"
  )

  lines <- c(lines, "", "matrix[N, J] mu;  // expected mixture isotope value")
  if (needs_proc) lines <- c(lines, "matrix[N, J] proc_var;  // source/TDF variance propagated into the mixture")

  lines <- c(lines, "for (i in 1:N) {", indent("for (j in 1:J) {  // mu[i, j] = proportion-weighted source + TDF means"))
  body <- character(0)
  if (spec$has_conc_dep) {
    body <- c(body,
      "real denom = 0;  // concentration-weighted normaliser",
      "for (k in 1:K) denom += p[i, k] * conc[k, j];",
      "real mij = 0;  // concentration-weighted expected isotope value"
    )
    if (needs_proc) body <- c(body, "real vij = 1e-8;  // concentration-weighted propagated variance")
    body <- c(body,
      "for (k in 1:K) {  // accumulate each source's concentration-weighted contribution",
      indent("real pk = p[i, k] * conc[k, j] / denom;  // concentration-adjusted proportion"),
      indent("mij += pk * (source_mean[k, j] + tdf_mean[k, j]);"),
      if (needs_proc) indent("vij += square(pk) * (square(source_sd[k, j]) + square(tdf_sd[k, j]));"),
      "}",
      "mu[i, j] = mij;"
    )
    if (needs_proc) body <- c(body, "proc_var[i, j] = vij;")
  } else {
    body <- c(body, "real mij = 0;  // expected isotope value")
    if (needs_proc) body <- c(body, "real vij = 1e-8;  // source/TDF variance propagated into the mixture")
    body <- c(body,
      "for (k in 1:K) {  // accumulate each source's contribution, weighted by its proportion",
      indent("mij += p[i, k] * (source_mean[k, j] + tdf_mean[k, j]);"),
      if (needs_proc) indent("vij += square(p[i, k]) * (square(source_sd[k, j]) + square(tdf_sd[k, j]));"),
      "}",
      "mu[i, j] = mij;"
    )
    if (needs_proc) body <- c(body, "proc_var[i, j] = vij;")
  }
  lines <- c(lines, indent(indent(body)), indent("}"), "}")

  if (needs_source_cor) {
    lines <- c(lines, "",
      "array[K] matrix[J, J] L_source_cov;  // per-source isotope covariance (Cholesky factor)",
      "array[K] matrix[J, J] source_cov;  // per-source isotope covariance (variance + correlation)",
      "for (k in 1:K) {",
      indent("L_source_cov[k] = diag_pre_multiply(to_vector(source_sd[k]), L_source_corr[k]);"),
      indent("source_cov[k] = multiply_lower_tri_self_transpose(L_source_cov[k]);"),
      "}"
    )
  }

  if (needs_mv) {
    cov_body <- if (spec$has_conc_dep) {
      c("real denom1 = 0;  // concentration-weighted normaliser for j1",
        "for (k in 1:K) denom1 += p[i, k] * conc[k, j1];",
        "real denom2 = 0;  // concentration-weighted normaliser for j2",
        "for (k in 1:K) denom2 += p[i, k] * conc[k, j2];",
        "real cij = 0;  // off-diagonal process covariance between isotopes j1, j2",
        "for (k in 1:K) {",
        indent("real pk1 = p[i, k] * conc[k, j1] / denom1;  // concentration-adjusted proportion for j1"),
        indent("real pk2 = p[i, k] * conc[k, j2] / denom2;  // concentration-adjusted proportion for j2"),
        indent("cij += pk1 * pk2 * source_cov[k][j1, j2];"),
        "}",
        "Omega[i][j1, j2] = cij;",
        "Omega[i][j2, j1] = cij;"
      )
    } else {
      c("real cij = 0;  // off-diagonal process covariance between isotopes j1, j2",
        "for (k in 1:K) {",
        indent("cij += square(p[i, k]) * source_cov[k][j1, j2];"),
        "}",
        "Omega[i][j1, j2] = cij;",
        "Omega[i][j2, j1] = cij;"
      )
    }
    lines <- c(lines, "",
      "array[N] matrix[J, J] Omega;  // full source/TDF covariance propagated into the mixture",
      "for (i in 1:N) {",
      indent("Omega[i] = diag_matrix(to_vector(proc_var[i]));  // diagonal: per-isotope process variance"),
      indent("for (j1 in 1:(J - 1)) {"),
      indent(indent("for (j2 in (j1 + 1):J) {")),
      indent(indent(indent(cov_body))),
      indent(indent("}")),
      indent("}"),
      "}"
    )

    lines <- c(lines, "",
      "array[N] matrix[J, J] L_Sigma;  // mixture covariance (Cholesky factor)",
      "for (i in 1:N) {",
      if (needs_resid_prop) {
        indent("L_Sigma[i] = diag_pre_multiply(sqrt(resid_prop), cholesky_decompose(Omega[i]));  // MixSIAR residual-error factor scales process covariance")
      } else {
        indent("L_Sigma[i] = cholesky_decompose(Omega[i]);")
      },
      "}"
    )
  }

  if (needs_resid_cor) {
    lines <- c(lines, "",
      "matrix[J, J] L_Sigma_resid = diag_pre_multiply(sigma, Lcorr_resid);  // shared residual covariance (Cholesky factor)"
    )
  }

  lines
}

## ---- model ---------------------------------------------------------------

#' Build the `model` block: a prior statement for every parameter (looked
#' up from `prior_df` via `select_prior()` -- `p_global`, `beta` if
#' `spec$P > 0`, each group-level term's `sd_re_*`/`Lcorr_re_*`, the
#' error-structure-specific parameters, and, for raw source/TDF data, their
#' estimated means/SDs/correlations), the raw-source and raw-TDF submodel
#' likelihoods (if applicable), and finally the mixture likelihood --
#' whose exact form (multivariate with full process+residual covariance,
#' multivariate with shared residual covariance, per-isotope normal with a
#' process-and/or-resid_prop-scaled SD, or plain per-isotope normal with
#' `sigma`) is chosen from `needs_mv`/`needs_resid_cor`/`needs_proc`/
#' `needs_resid_prop`.
#'
#' @param spec A `bsimms_spec` (see `build_bsimms_spec()`).
#' @param prior_df A merged `bsimms_prior` data frame (as returned by
#'   `merge_bsimms_prior()`).
#' @param needs_sigma,needs_resid_prop,needs_proc,needs_source_cor,needs_mv,needs_resid_cor
#'   See `bsimms_stancode_from_spec()`.
#' @return Character vector of Stan code lines.
#' @noRd
stan_model_lines <- function(spec, prior_df, needs_sigma, needs_resid_prop, needs_proc, needs_source_cor, needs_mv, needs_resid_cor) {
  lines <- character(0)

  lines <- c(lines, "// prior: global (population-average) source proportions",
             "p_global ~ dirichlet(alpha_dirichlet);")

  if (spec$P > 0) {
    lines <- c(lines, "", "// priors: fixed-effect covariate slopes")
    for (p in seq_len(spec$P)) {
      pr <- select_prior(prior_df, "b", coef = spec$fixed_names[p])
      lines <- c(lines, sprintf("beta[%d] ~ %s;  // %s", p, pr, spec$fixed_names[p]))
    }
  }

  if (length(spec$re_terms) > 0) {
    lines <- c(lines, "", "// priors: group-level effects")
    for (re in spec$re_terms) {
      size <- length(re$term_names) * spec$D
      lines <- c(lines, sprintf("to_vector(z_re_%s) ~ std_normal();", re$label))
      lines <- c(lines, sprintf("sd_re_%s ~ %s;", re$label, select_prior(prior_df, "sd", group = re$group)))
      if (size > 1) {
        lines <- c(lines, sprintf("Lcorr_re_%s ~ %s;", re$label, select_prior(prior_df, "cor", group = re$group)))
      }
    }
  }

  if (needs_sigma) {
    lines <- c(lines, "", "// priors: error term")
    for (j in seq_along(spec$isotope_names)) {
      lines <- c(lines, sprintf("sigma[%d] ~ %s;  // %s", j, select_prior(prior_df, "sigma", resp = spec$isotope_names[j]), spec$isotope_names[j]))
    }
  }

  if (needs_resid_cor) {
    lines <- c(lines, "", "// prior: residual-error correlation (Cholesky factor)",
               sprintf("Lcorr_resid ~ %s;", select_prior(prior_df, "resid_cor")))
  }

  if (needs_resid_prop) {
    lines <- c(lines, "", "// priors: MixSIAR residual-error factor (scales process variance)")
    for (j in seq_along(spec$isotope_names)) {
      lines <- c(lines, sprintf("resid_prop[%d] ~ %s;  // %s", j, select_prior(prior_df, "resid_prop", resp = spec$isotope_names[j]), spec$isotope_names[j]))
    }
  }

  if (spec$source$mode == "raw") {
    lines <- c(lines, "", "// priors + sub-model: raw source data (source- and isotope-specific)")
    for (k in seq_along(spec$source_names)) {
      for (j in seq_along(spec$isotope_names)) {
        pr_mean <- select_prior(prior_df, "source_mean", resp = spec$isotope_names[j], group = spec$source_names[k])
        pr_sd   <- select_prior(prior_df, "source_sd",   resp = spec$isotope_names[j], group = spec$source_names[k])
        lines <- c(lines, sprintf("source_mean[%d, %d] ~ %s;  // %s, %s", k, j, pr_mean, spec$source_names[k], spec$isotope_names[j]))
        lines <- c(lines, sprintf("source_sd[%d, %d] ~ %s;", k, j, pr_sd))
      }
    }

    if (needs_source_cor) {
      lines <- c(lines, "", "// priors: per-source isotope correlation (Cholesky factor)")
      for (k in seq_along(spec$source_names)) {
        pr_cor <- select_prior(prior_df, "source_cor", group = spec$source_names[k])
        lines <- c(lines, sprintf("L_source_corr[%d] ~ %s;  // %s", k, pr_cor, spec$source_names[k]))
      }
      lines <- c(lines,
        "",
        "// likelihood for each raw source measurement (isotopes jointly, per-source correlation)",
        "for (n in 1:N_source_raw) {",
        indent("to_vector(source_raw[n]) ~ multi_normal_cholesky(to_vector(source_mean[source_idx[n]]), L_source_cov[source_idx[n]]);"),
        "}"
      )
    } else {
      lines <- c(lines,
        "for (n in 1:N_source_raw) {  // likelihood for each raw source measurement",
        indent("source_raw[n] ~ normal(source_mean[source_idx[n]], source_sd[source_idx[n]]);"),
        "}"
      )
    }
  }

  if (spec$tdf$mode == "raw") {
    lines <- c(lines, "", "// priors + sub-model: raw TDF data (source- and isotope-specific)")
    for (k in seq_along(spec$source_names)) {
      for (j in seq_along(spec$isotope_names)) {
        pr_mean <- select_prior(prior_df, "tdf_mean", resp = spec$isotope_names[j], group = spec$source_names[k])
        pr_sd   <- select_prior(prior_df, "tdf_sd",   resp = spec$isotope_names[j], group = spec$source_names[k])
        lines <- c(lines, sprintf("tdf_mean[%d, %d] ~ %s;  // %s, %s", k, j, pr_mean, spec$source_names[k], spec$isotope_names[j]))
        lines <- c(lines, sprintf("tdf_sd[%d, %d] ~ %s;", k, j, pr_sd))
      }
    }
    lines <- c(lines,
      "for (n in 1:N_tdf_raw) {  // likelihood for each raw TDF measurement",
      indent("tdf_raw[n] ~ normal(tdf_mean[tdf_idx[n]], tdf_sd[tdf_idx[n]]);"),
      "}"
    )
  }

  lines <- c(lines, "", "// mixture likelihood")
  lik <- if (needs_mv) {
    c("for (i in 1:N) {  // isotopes jointly, full process (+ residual) covariance",
      indent("to_vector(y[i]) ~ multi_normal_cholesky(to_vector(mu[i]), L_Sigma[i]);"),
      "}")
  } else if (needs_resid_cor) {
    c("for (i in 1:N) {  // isotopes jointly, shared residual covariance",
      indent("to_vector(y[i]) ~ multi_normal_cholesky(to_vector(mu[i]), L_Sigma_resid);"),
      "}")
  } else if (needs_proc && needs_resid_prop) {
    c("for (i in 1:N) {",
      indent("for (j in 1:J) {"),
      indent(indent("y[i, j] ~ normal(mu[i, j], sqrt(proc_var[i, j] * resid_prop[j]));")),
      indent("}"),
      "}")
  } else if (needs_proc) {
    c("for (i in 1:N) {",
      indent("for (j in 1:J) {"),
      indent(indent("y[i, j] ~ normal(mu[i, j], sqrt(proc_var[i, j]));")),
      indent("}"),
      "}")
  } else {
    c("for (i in 1:N) {",
      indent("to_vector(y[i]) ~ normal(to_vector(mu[i]), sigma);"),
      "}")
  }
  c(lines, lik)
}

## ---- generated quantities -------------------------------------------------

#' Build the `generated quantities` block: `log_lik` (per-mixture-sample log
#' density, for LOO/WAIC model comparison) and `y_rep` (posterior
#' predictive draws), matching whichever mixture-likelihood form
#' `stan_model_lines()` used (selected the same way, from `needs_mv`/
#' `needs_resid_cor`/`needs_proc`/`needs_resid_prop`): a multivariate
#' `*_lpdf`/`*_rng` pair when `needs_mv` or `needs_resid_cor`, otherwise a
#' per-isotope loop accumulating `log_lik` and drawing `y_rep` from
#' `normal_lpdf`/`normal_rng` with the SD implied by `needs_proc`/
#' `needs_resid_prop` (process-and/or-resid_prop-scaled, or plain `sigma`).
#'
#' @param spec A `bsimms_spec` (see `build_bsimms_spec()`).
#' @param needs_sigma,needs_resid_prop,needs_proc,needs_mv,needs_resid_cor
#'   See `bsimms_stancode_from_spec()`.
#' @return Character vector of Stan code lines.
#' @noRd
stan_generated_quantities_lines <- function(spec, needs_sigma, needs_resid_prop, needs_proc, needs_mv, needs_resid_cor) {
  if (needs_mv) {
    return(c(
      "vector[N] log_lik;  // joint log density per mixture sample (for loo)",
      "matrix[N, J] y_rep;",
      "for (i in 1:N) {",
      indent("log_lik[i] = multi_normal_cholesky_lpdf(to_vector(y[i]) | to_vector(mu[i]), L_Sigma[i]);"),
      indent("y_rep[i] = to_row_vector(multi_normal_cholesky_rng(to_vector(mu[i]), L_Sigma[i]));"),
      "}"
    ))
  }

  if (needs_resid_cor) {
    return(c(
      "vector[N] log_lik;  // joint log density per mixture sample (for loo)",
      "matrix[N, J] y_rep;",
      "for (i in 1:N) {",
      indent("log_lik[i] = multi_normal_cholesky_lpdf(to_vector(y[i]) | to_vector(mu[i]), L_Sigma_resid);"),
      indent("y_rep[i] = to_row_vector(multi_normal_cholesky_rng(to_vector(mu[i]), L_Sigma_resid));"),
      "}"
    ))
  }

  sd_expr <- if (needs_proc && needs_resid_prop) {
    "sqrt(proc_var[i, j] * resid_prop[j])"
  } else if (needs_proc) {
    "sqrt(proc_var[i, j])"
  } else {
    "sigma[j]"
  }
  c(
    "vector[N] log_lik;  // joint log density per mixture sample (isotopes summed; for loo)",
    "matrix[N, J] y_rep;",
    "for (i in 1:N) {",
    indent("real ll_i = 0;"),
    indent("for (j in 1:J) {"),
    indent(indent(sprintf("real sd_ij = %s;", sd_expr))),
    indent(indent("ll_i += normal_lpdf(y[i, j] | mu[i, j], sd_ij);")),
    indent(indent("y_rep[i, j] = normal_rng(mu[i, j], sd_ij);")),
    indent("}"),
    indent("log_lik[i] = ll_i;"),
    "}"
  )
}
