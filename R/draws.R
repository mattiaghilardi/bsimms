#' Extract posterior draws from a `bsimms` fit
#'
#' Backend-agnostic accessor returning a `posterior::draws_array` regardless
#' of whether the model was fit with `cmdstanr` or `rstan`.
#'
#' @param object A `bsimms_fit` object (as returned by [bsimm()]).
#' @param variable Optional character vector of variable names (or
#'   `posterior`-style selectors) to extract; `NULL` extracts everything.
#' @return A `posterior::draws_array`.
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
#' bsimms_draws(fit, variable = "p_global")
#' }
bsimms_draws <- function(object, variable = NULL) {
  if (!inherits(object, "bsimms_fit")) {
    cli::cli_abort("{.arg object} must be a {.cls bsimms_fit} object.", call = NULL)
  }
  if (object$backend == "cmdstanr") {
    object$fit$draws(variables = variable)
  } else {
    d <- posterior::as_draws_array(object$fit)
    if (!is.null(variable)) d <- posterior::subset_draws(d, variable = variable)
    d
  }
}

#' Extract posterior draws from a `bsimms` fit as a `posterior::draws_matrix`
#' (chains collapsed into the draws dimension), the 2-D form most internal
#' post-processing helpers operate on (as opposed to `bsimms_draws()`'s 3-D
#' `draws_array`).
#'
#' @param object A `bsimms_fit` object.
#' @param variable Optional character vector of variable names (or
#'   `posterior`-style selectors) to extract; `NULL` extracts everything.
#' @return A `posterior::draws_matrix`.
#' @noRd
draws_matrix <- function(object, variable = NULL) {
  posterior::as_draws_matrix(bsimms_draws(object, variable = variable))
}

#' Randomly subset a `posterior::draws_matrix`'s rows (draws) down to
#' `ndraws`.
#'
#' @param dm A `posterior::draws_matrix` (as returned by `draws_matrix()`).
#' @param ndraws Number of draws to keep, or `NULL` (default) to keep all.
#' @return `dm`, subset to `ndraws` randomly chosen rows, unchanged if
#'   `ndraws` is `NULL`.
#' @noRd
subset_ndraws <- function(dm, ndraws) {
  if (is.null(ndraws)) {
    return(dm)
  }
  n_draws <- nrow(dm)
  if (!is.numeric(ndraws) || length(ndraws) != 1 || ndraws < 1 || ndraws != round(ndraws)) {
    cli::cli_abort("{.arg ndraws} must be a single positive integer.", call = NULL)
  }
  if (ndraws > n_draws) {
    cli::cli_abort(
      "{.arg ndraws} ({ndraws}) cannot exceed the number of posterior draws available ({n_draws}).",
      call = NULL
    )
  }
  dm[sample.int(n_draws, ndraws), , drop = FALSE]
}

#' Extract a rectangular `(n_draws x dim1)` or `(n_draws x dim1 x dim2)`
#' array of draws for a Stan parameter named e.g. `beta[p,d]` or `sigma[j]`
#' (i.e. reshapes the matching flat `prefix[i]`/`prefix[i,j]` columns of a
#' draws matrix into an indexed array).
#'
#' @param dm A `posterior::draws_matrix` (as returned by `draws_matrix()`).
#' @param prefix Stan parameter name, e.g. `"beta"` or `"sigma"`.
#' @param dim1 Size of the parameter's first index.
#' @param dim2 Size of the parameter's second index, if any (`NULL`,
#'   default, for a 1-D/vector parameter).
#' @return A numeric array: `[n_draws, dim1]` if `dim2` is `NULL`, otherwise
#'   `[n_draws, dim1, dim2]`.
#' @noRd
extract_array_draws <- function(dm, prefix, dim1, dim2 = NULL) {
  n_draws <- nrow(dm)
  if (is.null(dim2)) {
    nm <- sprintf("%s[%d]", prefix, seq_len(dim1))
    missing_nm <- setdiff(nm, colnames(dm))
    if (length(missing_nm) > 0) {
      cli::cli_abort("Could not find draws for: {.field {missing_nm}}.", call = NULL)
    }
    matrix(as.numeric(dm[, nm, drop = FALSE]), nrow = n_draws)
  } else {
    arr <- array(NA_real_, dim = c(n_draws, dim1, dim2))
    for (a in seq_len(dim1)) {
      nm <- sprintf("%s[%d,%d]", prefix, a, seq_len(dim2))
      missing_nm <- setdiff(nm, colnames(dm))
      if (length(missing_nm) > 0) {
        cli::cli_abort("Could not find draws for: {.field {missing_nm}}.", call = NULL)
      }
      arr[, a, ] <- dm[, nm, drop = FALSE]
    }
    arr
  }
}

#' Extract a rectangular `[n_draws, dim1, dim2, dim2]` array of draws for a
#' Stan parameter that is an array of square matrices, e.g.
#' `array[K] cholesky_factor_corr[J] L_source_corr`, flattened by Stan into
#' `prefix[k,j1,j2]`-named columns (i.e. `extract_array_draws()`'s 1-/2-index
#' extraction generalised to this 3-index case).
#'
#' @param dm A `posterior::draws_matrix` (as returned by `draws_matrix()`).
#' @param prefix Stan parameter name, e.g. `"L_source_corr"`.
#' @param dim1 Size of the parameter's array index (e.g. `K`).
#' @param dim2 Size of each matrix's row/column index (e.g. `J`).
#' @return A numeric array `[n_draws, dim1, dim2, dim2]`.
#' @noRd
extract_array_of_matrices <- function(dm, prefix, dim1, dim2) {
  n_draws <- nrow(dm)
  arr <- array(NA_real_, dim = c(n_draws, dim1, dim2, dim2))
  for (a in seq_len(dim1)) {
    for (b in seq_len(dim2)) {
      nm <- sprintf("%s[%d,%d,%d]", prefix, a, b, seq_len(dim2))
      missing_nm <- setdiff(nm, colnames(dm))
      if (length(missing_nm) > 0) {
        cli::cli_abort("Could not find draws for: {.field {missing_nm}}.", call = NULL)
      }
      arr[, a, b, ] <- dm[, nm, drop = FALSE]
    }
  }
  arr
}

#' Reshape a draws array into long format
#'
#' Reshapes a `[n_draws, n_obs, n_var]` array (as returned by
#' [posterior_proportions()], [posterior_epred.bsimms_fit()], or
#' [posterior_predict.bsimms_fit()]), with variable names attached as the
#' 3rd dimension's `dimnames`, into a long-format data frame with one row
#' per (draw, observation, variable) triple, ready for custom plots or
#' summaries (e.g. with `ggplot2`/`ggdist`).
#'
#' @param arr A numeric `[n_draws, n_obs, n_var]` array, with variable
#'   names attached as the 3rd dimension's `dimnames`.
#' @param var_col Name to give the variable-identity column (default
#'   `"variable"`), e.g. `"source"` for a [posterior_proportions()] array
#'   or `"isotope"` for a [posterior_epred.bsimms_fit()] array.
#' @param value_col Name to give the value column (default `"value"`).
#' @return A long-format data frame with columns `draw` (posterior draw
#'   index), `row` (observation index, into `newdata` or the fitted
#'   mixture samples), `<var_col>`, and `<value_col>`.
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
#' draws_long(p_arr, var_col = "source", value_col = "proportion")
#' }
draws_long <- function(arr, var_col = "variable", value_col = "value") {
  if (length(dim(arr)) != 3 || is.null(dimnames(arr)[[3]])) {
    cli::cli_abort(
      "{.arg arr} must be a `[n_draws, n_obs, n_var]` array with variable names attached as the 3rd dimension's dimnames.",
      call = NULL
    )
  }
  if (identical(var_col, value_col)) {
    cli::cli_abort("{.arg var_col} and {.arg value_col} must be different.", call = NULL)
  }
  n_draws <- dim(arr)[1]
  n_obs <- dim(arr)[2]
  var_names <- dimnames(arr)[[3]]

  out <- data.frame(
    draw = rep(seq_len(n_draws), times = n_obs * length(var_names)),
    row = rep(rep(seq_len(n_obs), each = n_draws), times = length(var_names))
  )
  out[[var_col]] <- rep(var_names, each = n_draws * n_obs)
  out[[value_col]] <- as.vector(arr)
  out
}
