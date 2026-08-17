#' Generate the Stan data list for a `bsimms` model
#'
#' Builds the named list of data passed to Stan, matching exactly the
#' `data` block produced by [make_stancode()] for the same arguments.
#'
#' @inheritParams make_stancode
#' @return A named list suitable for `cmdstanr::sample(data = ...)` or
#'   `rstan::sampling(data = ...)`. Class `bsimms_standata`.
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
#' sdata <- make_standata(
#'   ~ 1 + (1 | Region), mixture_data = mixture_data, source_data = source_data,
#'   tdf_data = tdf_data, isotope_names = c("d13C", "d15N")
#' )
#' str(sdata, max.level = 1)
make_standata <- function(formula, mixture_data, source_data, tdf_data, isotope_names,
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
  out <- standata_from_spec(spec, prior_df)
  class(out) <- c("bsimms_standata", "list")
  out
}

#' Build the named list of Stan data for a model, matching exactly the
#' `data` block declarations produced by `stan_data_lines()` for the same
#' `spec`: dimensions/design data, `alpha_dirichlet` (looked up from
#' `prior_df` via `select_prior()`, one `"p_global"` concentration per
#' source), `P`/`X` (only if `spec$P > 0`), one group of entries per
#' group-level term in `spec$re_terms`, source and TDF data (raw replicate
#' samples or means/SDs, depending on `spec$source`/`spec$tdf`'s `mode`),
#' and, if `spec$has_conc_dep`, `conc`.
#'
#' @param spec A `bsimms_spec` (see `build_bsimms_spec()`).
#' @param prior_df A merged `bsimms_prior` data frame (as returned by
#'   `merge_bsimms_prior()`).
#' @return A named list of Stan data, unclassed (`make_standata()` adds the
#'   `bsimms_standata` class).
#' @noRd
standata_from_spec <- function(spec, prior_df) {
  alpha_dirichlet <- vapply(
    spec$source_names,
    function(s) as.numeric(select_prior(prior_df, "p_global", group = s)),
    numeric(1)
  )
  out <- list(
    N = spec$N, J = spec$J, K = spec$K, D = spec$D,
    V = spec$V, y = spec$y,
    alpha_dirichlet = unname(alpha_dirichlet)
  )

  if (spec$P > 0) {
    out$P <- spec$P
    out$X <- spec$X
  }

  for (re in spec$re_terms) {
    out[[paste0("N_re_", re$label)]] <- length(re$group_levels)
    out[[paste0("M_re_", re$label)]] <- length(re$term_names)
    out[[paste0("Z_re_", re$label)]] <- re$Z
    out[[paste0("grp_re_", re$label)]] <- re$group_idx
  }

  if (spec$source$mode == "raw") {
    out$N_source_raw <- spec$source$n
    out$source_idx <- spec$source$source_idx
    out$source_raw <- spec$source$Y
  } else {
    out$source_mean_data <- spec$source$mean
    out$source_sd_data <- spec$source$sd
  }

  if (spec$tdf$mode == "raw") {
    out$N_tdf_raw <- spec$tdf$n
    out$tdf_idx <- spec$tdf$source_idx
    out$tdf_raw <- spec$tdf$Y
  } else {
    out$tdf_mean_data <- spec$tdf$mean
    out$tdf_sd_data <- spec$tdf$sd
  }

  if (spec$has_conc_dep) {
    out$conc <- spec$conc
  }

  out
}
