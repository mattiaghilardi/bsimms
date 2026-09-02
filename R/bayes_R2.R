#' Bayesian R-squared
#'
#' Computes a Bayesian R-squared (Gelman, Goodrich, Gabry, and Vehtari 2019)
#' for one or more isotopes. If `object$criteria$bayes_R2` was already
#' cached via [add_criterion()], it is reused (subset to `resp`) instead of
#' recomputed.
#'
#' @param object A `bsimms_fit` object (as returned by [bsimm()]).
#' @param resp Optional character vector of isotope names (a subset of
#'   `isotope_names`) to compute R-squared for. `NULL` (default) computes it for
#'   every isotope.
#' @param summary Logical; return a summary data frame (default) or the raw
#'   posterior draws.
#' @param robust Logical; if `FALSE` (default) summarise the central
#'   tendency/spread with `mean`/`sd`, if `TRUE` use `median`/`mad` instead.
#' @param probs Quantiles to include in the summary (default 2.5%/97.5%).
#' @param ... Currently unused.
#' @return If `summary = TRUE`, a data frame with one row per isotope
#'   (`isotope`, and one column per summary measure, named as in
#'   [posterior::summarise_draws()]: `mean`, `sd`, `q2.5`, `q97.5`, or with
#'   `robust = TRUE`, `median`, `mad`, `q2.5`, `q97.5`). If `FALSE`, an
#'   `n_draws x length(resp)` matrix of R-squared draws (column names = isotope
#'   names).
#' @references Gelman, A., Goodrich, B., Gabry, J., & Vehtari, A. (2019).
#'   R-squared for Bayesian regression models. *The American Statistician*,
#'   73(3), 307-309. \doi{10.1080/00031305.2018.1549100}.
#'   ([Preprint](https://acris.aalto.fi/ws/portalfiles/portal/34206843/bayes_R2_v3.pdf))
#' @exportS3Method rstantools::bayes_R2
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
#' rstantools::bayes_R2(fit)
#' }
bayes_R2.bsimms_fit <- function(object, resp = NULL, summary = TRUE, robust = FALSE,
                                 probs = c(0.025, 0.975), ...) {
  spec <- object$spec
  resp <- if (is.null(resp)) {
    spec$isotope_names
  } else {
    rlang::arg_match(resp, values = spec$isotope_names, multiple = TRUE)
  }

  cached <- object$criteria$bayes_R2
  if (!is.null(cached)) {
    r2 <- cached[, resp, drop = FALSE]
  } else {
    dm <- draws_matrix(object, variable = "mu")
    mu_arr <- extract_array_draws(dm, "mu", spec$N, spec$J) # n_draws x N x J
    y <- object$standata$y

    j_idx <- match(resp, spec$isotope_names)
    r2 <- matrix(NA_real_, nrow(mu_arr), length(resp))
    colnames(r2) <- resp
    for (k in seq_along(j_idx)) {
      fit_j <- mu_arr[, , j_idx[k]] # n_draws x N
      resid_j <- sweep(-fit_j, 2, y[, j_idx[k]], "+") # y[i] - fit[s, i]
      var_fit <- apply(fit_j, 1, stats::var)
      var_res <- apply(resid_j, 1, stats::var)
      r2[, k] <- var_fit / (var_fit + var_res)
    }
  }

  if (!summary) return(r2)
  measures <- c(
    if (robust) list("median", "mad") else list("mean", "sd"),
    list(~ posterior::quantile2(.x, probs = probs))
  )
  out <- as.data.frame(do.call(
    posterior::summarise_draws, c(list(posterior::as_draws_matrix(r2)), measures)
  ))
  names(out)[names(out) == "variable"] <- "isotope"
  rownames(out) <- NULL
  out
}
