mixture_data <- data.frame(
  d13C = c(-20, -21, -19, -22, -20.5, -21.5),
  d15N = c(10, 11, 9, 12, 10.5, 11.5),
  elevation = c(100, 120, 110, 300, 320, 310),
  Region = factor(rep(c("A", "B"), each = 3))
)
source_data <- data.frame(
  Source = c("Beaver", "Deer", "Otter"),
  d13C_mean = c(-26, -14, -20), d13C_sd = c(1, 1, 1),
  d15N_mean = c(4, 6, 3), d15N_sd = c(1, 1, 1)
)
tdf_data <- data.frame(
  Source = c("Beaver", "Deer", "Otter"),
  d13C_mean = c(0.8, 0.9, 1.0), d13C_sd = c(0.2, 0.2, 0.2),
  d15N_mean = c(3.2, 3.3, 3.4), d15N_sd = c(0.3, 0.3, 0.3)
)

skip_if_not_installed("cmdstanr")
skip_if_not_installed("rstantools")

fit <- bsimm(
  formula = ~ elevation + (1 | Region), mixture_data = mixture_data, source_data = source_data,
  tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), source_means_sds = TRUE,
  chains = 1, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0,
  show_messages = FALSE, show_exceptions = FALSE
)

raw_source_data <- data.frame(
  Source = rep(c("Beaver", "Deer"), each = 4),
  d13C = c(-25.5, -26.2, -24.8, -25.9, -18.1, -17.5, -18.6, -17.9),
  d15N = c(4.8, 5.1, 4.6, 5.0, 7.9, 8.2, 7.6, 8.0)
)
raw_tdf_data <- tdf_data[tdf_data$Source %in% c("Beaver", "Deer"), ]
raw_fit <- bsimm(
  formula = ~ elevation + (1 | Region), mixture_data = mixture_data, source_data = raw_source_data,
  tdf_data = raw_tdf_data, isotope_names = c("d13C", "d15N"), source_means_sds = FALSE,
  chains = 1, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0,
  show_messages = FALSE, show_exceptions = FALSE
)

# posterior_proportions / fitted_proportions ---------------------------------

test_that("posterior_proportions returns an [n_draws, N, K] array summing to 1", {
  p_arr <- posterior_proportions(fit)
  expect_equal(dim(p_arr), c(100, 6, 3))
  expect_true(all(abs(apply(p_arr, c(1, 2), sum) - 1) < 1e-6))
})

test_that("fitted_proportions summarises proportions by source", {
  s <- fitted_proportions(fit)
  expect_s3_class(s, "data.frame")
  expect_equal(unique(s$source), c("Beaver", "Deer", "Otter"))
  expect_equal(nrow(s), 6 * 3)
  expect_true(all(c("mean", "sd", "q2.5", "q97.5") %in% names(s)))
})

test_that("fitted_proportions with summary = FALSE matches posterior_proportions", {
  expect_equal(fitted_proportions(fit, summary = FALSE), posterior_proportions(fit))
})

test_that("fitted_proportions's robust argument switches mean/sd for median/mad", {
  s <- fitted_proportions(fit, robust = TRUE)
  expect_true(all(c("median", "mad", "q2.5", "q97.5") %in% names(s)))
})

test_that("posterior_proportions's ndraws subsets the number of draws returned", {
  p_arr <- posterior_proportions(fit, ndraws = 10)
  expect_equal(dim(p_arr), c(10, 6, 3))
})

test_that("fitted_proportions's ndraws is forwarded to posterior_proportions", {
  s <- fitted_proportions(fit, ndraws = 10, summary = FALSE)
  expect_equal(dim(s), c(10, 6, 3))
})

test_that("posterior_proportions errors when ndraws exceeds the number of draws available", {
  expect_snapshot(error = TRUE, posterior_proportions(fit, ndraws = 1000))
})

test_that("posterior_proportions predicts at new covariate combinations", {
  newdata <- data.frame(elevation = c(100, 200), Region = factor(c("A", "B")))
  p_arr <- posterior_proportions(fit, newdata = newdata)
  expect_equal(dim(p_arr), c(100, 2, 3))
  expect_true(all(abs(apply(p_arr, c(1, 2), sum) - 1) < 1e-6))
})

test_that("posterior_proportions's ndraws also applies to new-data prediction", {
  newdata <- data.frame(elevation = c(100, 200), Region = factor(c("A", "B")))
  p_arr <- posterior_proportions(fit, newdata = newdata, ndraws = 10)
  expect_equal(dim(p_arr), c(10, 2, 3))
})

test_that("posterior_proportions errors by default (re_formula = NULL) when newdata is missing a group-level column", {
  newdata <- data.frame(elevation = c(100, 200))
  expect_snapshot(error = TRUE, posterior_proportions(fit, newdata = newdata))
})

test_that("posterior_proportions predicts at the population-average level via re_formula = NA", {
  newdata <- data.frame(elevation = c(100, 200))
  p_arr <- posterior_proportions(fit, newdata = newdata, re_formula = NA)
  expect_equal(dim(p_arr), c(100, 2, 3))
})

test_that("posterior_proportions predicts at the population-average level via re_formula = ~0", {
  newdata <- data.frame(elevation = c(100, 200))
  p_arr <- posterior_proportions(fit, newdata = newdata, re_formula = ~0)
  expect_equal(dim(p_arr), c(100, 2, 3))
})

test_that("posterior_proportions errors on unseen group levels in newdata", {
  newdata <- data.frame(elevation = 100, Region = factor("C"))
  expect_snapshot(error = TRUE, posterior_proportions(fit, newdata = newdata))
})

test_that("posterior_proportions treats NA in a grouping column as a new level", {
  newdata <- data.frame(elevation = 100, Region = factor(NA, levels = c("A", "B")))
  expect_snapshot(error = TRUE, posterior_proportions(fit, newdata = newdata))
})

test_that("posterior_proportions errors when newdata is missing a fixed-effect column", {
  newdata <- data.frame(Region = factor("A"))
  expect_snapshot(error = TRUE, posterior_proportions(fit, newdata = newdata))
})

test_that("posterior_proportions errors when re_formula names a term not in the model", {
  newdata <- data.frame(elevation = 100, Region = factor("A"))
  expect_snapshot(error = TRUE, posterior_proportions(fit, newdata = newdata, re_formula = ~ (1 | Foo)))
})

test_that("posterior_proportions errors when re_formula is not NULL/NA/a formula", {
  newdata <- data.frame(elevation = 100, Region = factor("A"))
  expect_snapshot(error = TRUE, posterior_proportions(fit, newdata = newdata, re_formula = "banana"))
})

# posterior_proportions with allow_new_levels / sample_new_levels -----------

test_that("posterior_proportions samples new levels when allow_new_levels = TRUE (uncertainty)", {
  newdata <- data.frame(elevation = 100, Region = factor("C", levels = c("A", "B", "C")))
  p_arr <- posterior_proportions(fit, newdata = newdata, allow_new_levels = TRUE, sample_new_levels = "uncertainty")
  expect_equal(dim(p_arr), c(100, 1, 3))
  expect_true(all(abs(apply(p_arr, c(1, 2), sum) - 1) < 1e-6))
})

test_that("posterior_proportions samples new levels when allow_new_levels = TRUE (gaussian)", {
  newdata <- data.frame(elevation = 100, Region = factor("C", levels = c("A", "B", "C")))
  p_arr <- posterior_proportions(fit, newdata = newdata, allow_new_levels = TRUE, sample_new_levels = "gaussian")
  expect_equal(dim(p_arr), c(100, 1, 3))
  expect_true(all(abs(apply(p_arr, c(1, 2), sum) - 1) < 1e-6))
})

test_that("rows sharing the same new, non-NA level get identical sampled draws", {
  newdata <- data.frame(
    elevation = c(100, 100),
    Region = factor(c("C", "C"), levels = c("A", "B", "C"))
  )
  p_arr <- posterior_proportions(fit, newdata = newdata, allow_new_levels = TRUE)
  expect_equal(p_arr[, 1, ], p_arr[, 2, ])
})

test_that("distinct NA rows get independently sampled draws", {
  newdata <- data.frame(
    elevation = c(100, 100),
    Region = factor(c(NA, NA), levels = c("A", "B"))
  )
  p_arr <- posterior_proportions(fit, newdata = newdata, allow_new_levels = TRUE, sample_new_levels = "gaussian")
  expect_false(isTRUE(all.equal(p_arr[, 1, ], p_arr[, 2, ])))
})

test_that("re_formula = NA ignores allow_new_levels entirely (term not conditioned on)", {
  newdata <- data.frame(elevation = 100, Region = factor("C", levels = c("A", "B", "C")))
  p_arr <- posterior_proportions(fit, newdata = newdata, re_formula = NA)
  expect_equal(dim(p_arr), c(100, 1, 3))
})

# posterior_proportions with re_formula and newdata = NULL -------------------

test_that("re_formula = NA gives population-average proportions for the fitted mixture data", {
  p_default <- posterior_proportions(fit)
  p_marginal <- posterior_proportions(fit, re_formula = NA)
  expect_equal(dim(p_marginal), dim(p_default))
  expect_false(isTRUE(all.equal(p_marginal, p_default)))
})

test_that("re_formula = ~0 matches re_formula = NA for the fitted mixture data", {
  expect_equal(
    posterior_proportions(fit, re_formula = ~0),
    posterior_proportions(fit, re_formula = NA)
  )
})

test_that("re_formula = NULL matches omitting re_formula for the fitted mixture data", {
  expect_equal(posterior_proportions(fit, re_formula = NULL), posterior_proportions(fit))
})

test_that("re_formula selecting all group-level terms matches the default for the fitted mixture data", {
  expect_equal(
    posterior_proportions(fit, re_formula = ~ (1 | Region)), posterior_proportions(fit),
    tolerance = 1e-6
  )
})

test_that("re_formula ignores allow_new_levels/sample_new_levels when newdata is NULL", {
  p1 <- posterior_proportions(fit, re_formula = NA, allow_new_levels = FALSE)
  p2 <- posterior_proportions(fit, re_formula = NA, allow_new_levels = TRUE, sample_new_levels = "gaussian")
  expect_equal(p1, p2)
})

# posterior_epred / fitted ----------------------------------------------------

test_that("posterior_epred returns an [n_draws, N, J] array for every isotope by default", {
  mu_arr <- rstantools::posterior_epred(fit)
  expect_equal(dim(mu_arr), c(100, 6, 2))
})

test_that("posterior_epred respects an explicit resp", {
  mu_arr <- rstantools::posterior_epred(fit, resp = "d15N")
  expect_equal(dim(mu_arr), c(100, 6, 1))
})

test_that("posterior_epred errors on an invalid resp", {
  expect_snapshot(error = TRUE, rstantools::posterior_epred(fit, resp = "banana"))
})

test_that("fitted summarises mu by isotope", {
  s <- fitted(fit)
  expect_s3_class(s, "data.frame")
  expect_equal(unique(s$isotope), c("d13C", "d15N"))
  expect_equal(nrow(s), 6 * 2)
})

test_that("fitted with summary = FALSE matches posterior_epred", {
  expect_equal(fitted(fit, summary = FALSE), rstantools::posterior_epred(fit))
})

test_that("posterior_epred's ndraws subsets the number of draws returned", {
  mu_arr <- rstantools::posterior_epred(fit, ndraws = 10)
  expect_equal(dim(mu_arr), c(10, 6, 2))
})

test_that("fitted's ndraws is forwarded to posterior_epred", {
  s <- fitted(fit, ndraws = 10, summary = FALSE)
  expect_equal(dim(s), c(10, 6, 2))
})

# posterior_predict / predict -------------------------------------------------

test_that("posterior_predict returns an [n_draws, N, J] array for every isotope by default", {
  yrep_arr <- rstantools::posterior_predict(fit)
  expect_equal(dim(yrep_arr), c(100, 6, 2))
})

test_that("posterior_predict respects an explicit resp", {
  yrep_arr <- rstantools::posterior_predict(fit, resp = "d13C")
  expect_equal(dim(yrep_arr), c(100, 6, 1))
})

test_that("predict summarises y_rep by isotope", {
  s <- predict(fit)
  expect_s3_class(s, "data.frame")
  expect_equal(unique(s$isotope), c("d13C", "d15N"))
  expect_equal(nrow(s), 6 * 2)
})

test_that("predict with summary = FALSE matches posterior_predict", {
  expect_equal(predict(fit, summary = FALSE), rstantools::posterior_predict(fit))
})

test_that("predict's probs argument controls the quantile columns", {
  s <- predict(fit, probs = c(0.1, 0.9))
  expect_true(all(c("q10", "q90") %in% names(s)))
})

test_that("posterior_predict's ndraws subsets the number of draws returned", {
  yrep_arr <- rstantools::posterior_predict(fit, ndraws = 10)
  expect_equal(dim(yrep_arr), c(10, 6, 2))
})

test_that("predict's ndraws is forwarded to posterior_predict", {
  s <- predict(fit, ndraws = 10, summary = FALSE)
  expect_equal(dim(s), c(10, 6, 2))
})

# posterior_epred / posterior_predict with re_formula and newdata = NULL ----

test_that("posterior_epred's re_formula gives population-average mu for the fitted mixture data", {
  mu_default <- rstantools::posterior_epred(fit)
  mu_marginal <- rstantools::posterior_epred(fit, re_formula = NA)
  expect_equal(dim(mu_marginal), dim(mu_default))
  expect_false(isTRUE(all.equal(mu_marginal, mu_default)))
})

test_that("posterior_predict's re_formula gives population-average y_rep for the fitted mixture data", {
  yrep_arr <- rstantools::posterior_predict(fit, re_formula = NA)
  expect_equal(dim(yrep_arr), c(100, 6, 2))
  expect_true(all(is.finite(yrep_arr)))
})

# posterior_epred / posterior_predict with newdata ---------------------------

test_that("posterior_epred reproduces fitted mu when newdata repeats the fitted mixture data", {
  mu_fitted <- rstantools::posterior_epred(fit)
  mu_new <- rstantools::posterior_epred(fit, newdata = mixture_data)
  expect_equal(mu_new, mu_fitted, tolerance = 1e-5)
})

test_that("posterior_epred reproduces fitted mu for raw source data too (mv-eligible model)", {
  mu_fitted <- rstantools::posterior_epred(raw_fit)
  mu_new <- rstantools::posterior_epred(raw_fit, newdata = mixture_data)
  expect_equal(mu_new, mu_fitted, tolerance = 1e-5)
})

test_that("posterior_epred's newdata prediction respects resp/re_formula/allow_new_levels/ndraws", {
  newdata <- data.frame(elevation = 100, Region = factor("C", levels = c("A", "B", "C")))
  mu_new <- rstantools::posterior_epred(
    fit,
    newdata = newdata, resp = "d15N", allow_new_levels = TRUE, ndraws = 10
  )
  expect_equal(dim(mu_new), c(10, 1, 1))
})

test_that("posterior_epred propagates posterior_proportions's newdata validation errors", {
  newdata <- data.frame(elevation = c(100, 200))
  expect_snapshot(error = TRUE, rstantools::posterior_epred(fit, newdata = newdata))
})

test_that("posterior_predict for new data adds observation noise (not just point predictions)", {
  newdata <- data.frame(elevation = c(100, 200), Region = factor(c("A", "B")))
  mu_new <- rstantools::posterior_epred(fit, newdata = newdata)
  yrep_new <- rstantools::posterior_predict(fit, newdata = newdata)
  expect_equal(dim(yrep_new), dim(mu_new))
  expect_true(all(is.finite(yrep_new)))
  expect_false(isTRUE(all.equal(yrep_new, mu_new)))
})

test_that("posterior_predict handles the multivariate (raw source, process + residual) branch for new data", {
  newdata <- data.frame(elevation = c(100, 200), Region = factor(c("A", "B")))
  yrep_new <- rstantools::posterior_predict(raw_fit, newdata = newdata)
  expect_equal(dim(yrep_new), c(100, 2, 2))
  expect_true(all(is.finite(yrep_new)))
})

test_that("fitted summarises mu for new data", {
  newdata <- data.frame(elevation = c(100, 200), Region = factor(c("A", "B")))
  s <- fitted(fit, newdata = newdata)
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), 2 * 2)
})

test_that("predict summarises y_rep for new data", {
  newdata <- data.frame(elevation = c(100, 200), Region = factor(c("A", "B")))
  s <- predict(fit, newdata = newdata)
  expect_s3_class(s, "data.frame")
  expect_equal(nrow(s), 2 * 2)
})

# posterior_proportions with a nested/crossed group-level term ---------------

nested_mixture_data <- data.frame(
  d13C = c(-20, -21, -19, -22, -20.5, -21.5, -20.2, -21.2),
  d15N = c(10, 11, 9, 12, 10.5, 11.5, 10.3, 11.3),
  Site = factor(rep(c("S1", "S2"), each = 4)),
  Individual = factor(c("I1", "I1", "I2", "I2", "I3", "I3", "I4", "I4"))
)
nested_source_data <- data.frame(
  Source = c("Beaver", "Deer"),
  d13C_mean = c(-25, -18), d13C_sd = c(1, 1),
  d15N_mean = c(5, 8), d15N_sd = c(1, 1)
)
nested_tdf_data <- data.frame(
  Source = c("Beaver", "Deer"),
  d13C_mean = c(1, 1.2), d13C_sd = c(0.2, 0.3),
  d15N_mean = c(3, 3.1), d15N_sd = c(0.4, 0.5)
)

nested_fit <- bsimm(
  formula = ~ 1 + (1 | Site / Individual), mixture_data = nested_mixture_data,
  source_data = nested_source_data, tdf_data = nested_tdf_data,
  isotope_names = c("d13C", "d15N"), source_means_sds = TRUE,
  chains = 1, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0,
  show_messages = FALSE, show_exceptions = FALSE
)

test_that("posterior_proportions correctly handles a nested group term with character newdata", {
  # row 1 of nested_mixture_data is Site = "S1", Individual = "I1"
  newdata_factor <- data.frame(Site = factor("S1"), Individual = factor("I1"))
  newdata_character <- data.frame(Site = "S1", Individual = "I1", stringsAsFactors = FALSE)
  p_fitted <- posterior_proportions(nested_fit)
  p_factor <- posterior_proportions(nested_fit, newdata = newdata_factor)
  p_character <- posterior_proportions(nested_fit, newdata = newdata_character)
  expect_equal(p_factor[, 1, ], p_fitted[, 1, ], tolerance = 1e-4)
  expect_equal(p_character[, 1, ], p_fitted[, 1, ], tolerance = 1e-4)
})

test_that("posterior_proportions errors by default when the inner nesting term's data is omitted", {
  expect_snapshot(
    error = TRUE,
    posterior_proportions(nested_fit, newdata = data.frame(Site = factor("S1")))
  )
})

test_that("re_formula selects the outer nesting term alone (site average over individuals)", {
  p <- posterior_proportions(nested_fit, newdata = data.frame(Site = factor("S1")), re_formula = ~ (1 | Site))
  expect_equal(dim(p), c(100, 1, 2))
})

test_that("re_formula = NA gives the pure population-average (neither nesting level)", {
  p <- posterior_proportions(nested_fit, newdata = data.frame(dummy = 1), re_formula = NA)
  expect_equal(dim(p), c(100, 1, 2))
})

test_that("posterior_proportions errors when the inner nesting level is supplied without the outer one", {
  expect_snapshot(
    error = TRUE,
    posterior_proportions(nested_fit, newdata = data.frame(Individual = factor("I1")))
  )
})

test_that("posterior_proportions errors on an unseen Individual/Site combination, referencing the actual columns", {
  # Individual "I3" only ever occurs with Site "S2" in nested_mixture_data
  expect_snapshot(
    error = TRUE,
    posterior_proportions(nested_fit, newdata = data.frame(Site = factor("S1"), Individual = factor("I3")))
  )
})

# posterior_proportions with a factor fixed effect ---------------------------

factor_fit <- bsimm(
  formula = ~Region, mixture_data = mixture_data, source_data = source_data,
  tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), source_means_sds = TRUE,
  chains = 1, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0,
  show_messages = FALSE, show_exceptions = FALSE
)

test_that("posterior_proportions handles a single-row newdata for a factor fixed effect", {
  # mixture_data row 1 is Region = "A"
  p_fitted <- posterior_proportions(factor_fit)
  p_factor <- posterior_proportions(factor_fit, newdata = data.frame(Region = factor("A", levels = c("A", "B"))))
  p_character <- posterior_proportions(factor_fit, newdata = data.frame(Region = "A"))
  p_underlevelled <- posterior_proportions(factor_fit, newdata = data.frame(Region = factor("A")))
  expect_equal(p_factor[, 1, ], p_fitted[, 1, ], tolerance = 1e-4)
  expect_equal(p_character[, 1, ], p_fitted[, 1, ], tolerance = 1e-4)
  expect_equal(p_underlevelled[, 1, ], p_fitted[, 1, ], tolerance = 1e-4)
})

test_that("posterior_proportions errors on a non-factor newdata column for a factor fixed effect", {
  expect_snapshot(error = TRUE, posterior_proportions(factor_fit, newdata = data.frame(Region = 1)))
})

test_that("posterior_proportions errors on an unseen factor level for a fixed effect", {
  expect_snapshot(
    error = TRUE,
    posterior_proportions(factor_fit, newdata = data.frame(Region = factor("C", levels = c("A", "B", "C"))))
  )
})
