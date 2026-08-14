#' Build the internal `bsimms_spec` model specification shared by
#' `make_stancode()`, `make_standata()`, `bsimms_get_prior()` and `bsimm()`:
#' parses and validates `formula`/`mixture_data`/`source_data`/`tdf_data`/
#' `conc_dep` (via `parse_bsimms_formula()`, `prep_mixture_isotopes()`,
#' `prep_iso_table()`, `prep_conc_dep()`), and derives every dimension
#' (`K`/`J`/`D`/`N`/`P`), matrix, and flag needed to generate Stan code/data
#' and the default priors, so those downstream steps never re-derive or
#' re-validate anything from the raw arguments themselves.
#'
#' @param formula,mixture_data,source_data,tdf_data,isotope_names,source_means_sds,tdf_means_sds,conc_dep,error_structure,source_col
#'   Same meaning as the corresponding arguments of `make_stancode()`.
#' @return A list of class `bsimms_spec` with elements: `formula`;
#'   `isotope_names`, `source_names` (sorted), `source_col`; dimensions `K`
#'   (sources), `J` (isotopes), `D` (`K - 1`, ILR dimensions), `N`
#'   (mixture samples), `P` (fixed-effect covariates, `(Intercept)`
#'   excluded); `V` (ILR basis, from `ilr_basis()`); `y` (mixture isotope
#'   matrix); `source` and `tdf` (as returned by `prep_iso_table()`); `conc`
#'   (as returned by `prep_conc_dep()`, or `NULL`) and `has_conc_dep`; `X`,
#'   `fixed_formula`, `fixed_names`, `fixed_frame` (fixed-effect design,
#'   from `parse_bsimms_formula()`, with `(Intercept)` dropped since the
#'   population-level baseline is `p_global`, not a `beta` intercept);
#'   `re_terms` (group-level terms, from `parse_bsimms_formula()`);
#'   `error_structure`; `source_means_sds`; `tdf_means_sds`.
#' @noRd
build_bsimms_spec <- function(formula,
                               mixture_data,
                               source_data,
                               tdf_data,
                               isotope_names,
                               source_means_sds = FALSE,
                               tdf_means_sds = TRUE,
                               conc_dep = FALSE,
                               error_structure = c("process_residual", "process_only", "residual_only"),
                               source_col = "Source") {
  error_structure <- rlang::arg_match(error_structure)

  if (!is.character(isotope_names) || length(isotope_names) < 1) {
    cli::cli_abort(
      "{.arg isotope_names} must be a non-empty character vector of mixture/source/TDF column names.",
      call = NULL
    )
  }

  source_names <- infer_source_names(source_data, tdf_data, source_col = source_col)
  K <- length(source_names)
  if (K < 2) cli::cli_abort("Need at least 2 sources.", call = NULL)
  J <- length(isotope_names)
  D <- K - 1L
  V <- ilr_basis(K)

  y <- prep_mixture_isotopes(mixture_data, isotope_names)
  N <- nrow(y)

  source <- prep_iso_table(source_data, isotope_names, source_names,
                            means_sds = source_means_sds, source_col = source_col, label = "source")
  tdf <- prep_iso_table(tdf_data, isotope_names, source_names,
                         means_sds = tdf_means_sds, source_col = source_col, label = "tdf")
  conc <- prep_conc_dep(source_data, isotope_names, source_names, source_means_sds,
                         conc_dep, source_col = source_col)

  pf <- parse_bsimms_formula(formula, mixture_data)
  if (nrow(pf$X) != N) {
    cli::cli_abort(
      c(
        "Design matrix built from {.arg formula} has {nrow(pf$X)} row{?s}, but the mixture isotope data has {N}.",
        "i" = "This is usually caused by missing values in the fixed-effect covariates: rows with {.code NA} are silently dropped when building the design matrix.",
        "i" = "Remove or impute the missing covariate values in {.arg mixture_data} before fitting."
      ),
      call = NULL
    )
  }
  for (re in pf$re_terms) {
    if (length(re$group_idx) != N) {
      cli::cli_abort(
        "Group-level term {.field {re$group}} has a different length than the mixture data.",
        call = NULL
      )
    }
  }

  # The population-level baseline is `p_global` (Dirichlet-distributed
  # simplex, see `default_bsimms_prior()`/`stancode.R`), not a `beta`
  # intercept, so drop `(Intercept)` from the fixed-effect design; any
  # remaining columns are pure covariate effects added on top of `p_global`
  # (in ILR space, via `ilr_global`). A formula with no other covariates
  # (e.g. `~ 1`) is a valid, intercept-only-via-`p_global` model (`P = 0`).
  keep <- pf$fixed_names != "(Intercept)"
  X <- pf$X[, keep, drop = FALSE]
  fixed_names <- pf$fixed_names[keep]

  structure(
    list(
      formula = formula,
      isotope_names = isotope_names,
      source_names = source_names,
      source_col = source_col,
      K = K, J = J, D = D, N = N, V = V,
      y = y,
      source = source,
      tdf = tdf,
      conc = conc,
      X = X,
      fixed_formula = pf$fixed_formula,
      fixed_names = fixed_names,
      fixed_frame = pf$fixed_frame,
      P = ncol(X),
      re_terms = pf$re_terms,
      error_structure = error_structure,
      source_means_sds = source_means_sds,
      tdf_means_sds = tdf_means_sds,
      has_conc_dep = !is.null(conc)
    ),
    class = "bsimms_spec"
  )
}
