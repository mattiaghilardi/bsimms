skip_if_no_backend <- function() {
  if (
    !requireNamespace("cmdstanr", quietly = TRUE) &&
      !requireNamespace("rstan", quietly = TRUE)
  ) {
    testthat::skip("cmdstanr or rstan required")
  }
}
