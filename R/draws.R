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
#' fit <- bsimm(
#'   ~ 1 + (1 | Region), mixture_data = mixture_data, source_data = source_data,
#'   tdf_data = tdf_data, isotope_names = c("d13C", "d15N"),
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
