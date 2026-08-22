mixture_data <- data.frame(
  d13C = c(-20, -21, -19, -22, -20.5, -21.5),
  d15N = c(10, 11, 9, 12, 10.5, 11.5),
  Region = factor(rep(c("A", "B"), each = 3))
)
source_data <- data.frame(
  Source = rep(c("Beaver", "Deer"), each = 3),
  d13C = c(-25, -24, -26, -18, -17, -19),
  d15N = c(5, 6, 4, 8, 9, 7)
)
tdf_data <- data.frame(
  Source = c("Beaver", "Deer"),
  d13C_mean = c(1, 1.2), d13C_sd = c(0.2, 0.3),
  d15N_mean = c(3, 3.1), d15N_sd = c(0.4, 0.5)
)

test_that("bsimms_draws errors on a non-bsimms_fit object", {
  expect_snapshot(error = TRUE, bsimms_draws(list()))
})

test_that("bsimms_draws returns a draws_array for a cmdstanr fit", {
  skip_if_not_installed("cmdstanr")
  fit <- bsimm(
    formula = ~1, mixture_data = mixture_data, source_data = source_data,
    tdf_data = tdf_data, isotope_names = c("d13C", "d15N"),
    chains = 1, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0
  )
  draws <- bsimms_draws(fit)
  expect_s3_class(draws, "draws_array")

  draws_p <- bsimms_draws(fit, variable = "p_global")
  expect_equal(posterior::variables(draws_p), c("p_global[1]", "p_global[2]"))
})

test_that("bsimms_draws returns a draws_array for an rstan fit", {
  skip_if_not_installed("rstan")
  fit <- suppressWarnings(bsimm(
    formula = ~1, mixture_data = mixture_data, source_data = source_data,
    tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), backend = "rstan",
    chains = 1, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0
  ))
  draws <- bsimms_draws(fit)
  expect_s3_class(draws, "draws_array")

  draws_p <- bsimms_draws(fit, variable = "p_global")
  expect_equal(posterior::variables(draws_p), c("p_global[1]", "p_global[2]"))
})

test_that("draws_matrix collapses chains into the draws dimension", {
  skip_if_not_installed("cmdstanr")
  fit <- bsimm(
    formula = ~1, mixture_data = mixture_data, source_data = source_data,
    tdf_data = tdf_data, isotope_names = c("d13C", "d15N"),
    chains = 2, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0
  )
  dm <- draws_matrix(fit, variable = "p_global")
  expect_s3_class(dm, "draws_matrix")
  expect_equal(nrow(dm), 200)
})

test_that("extract_array_draws reshapes a vector parameter", {
  skip_if_not_installed("cmdstanr")
  fit <- bsimm(
    formula = ~1, mixture_data = mixture_data, source_data = source_data,
    tdf_data = tdf_data, isotope_names = c("d13C", "d15N"),
    chains = 1, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0
  )
  dm <- draws_matrix(fit, variable = "p_global")
  arr <- extract_array_draws(dm, "p_global", dim1 = 2)
  expect_equal(dim(arr), c(100, 2))
  expect_equal(arr[, 1], as.numeric(dm[, "p_global[1]"]))
})

test_that("extract_array_draws reshapes a matrix parameter", {
  skip_if_not_installed("cmdstanr")
  fit <- bsimm(
    formula = ~ 1 + (1 | Region), mixture_data = mixture_data, source_data = source_data,
    tdf_data = tdf_data, isotope_names = c("d13C", "d15N"),
    chains = 1, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0
  )
  dm <- draws_matrix(fit, variable = "source_mean")
  arr <- extract_array_draws(dm, "source_mean", dim1 = 2, dim2 = 2)
  expect_equal(dim(arr), c(100, 2, 2))
  expect_equal(arr[, 1, 1], as.numeric(dm[, "source_mean[1,1]"]))
})

test_that("extract_array_draws errors when the requested draws are missing", {
  skip_if_not_installed("cmdstanr")
  fit <- bsimm(
    formula = ~1, mixture_data = mixture_data, source_data = source_data,
    tdf_data = tdf_data, isotope_names = c("d13C", "d15N"),
    chains = 1, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0
  )
  dm <- draws_matrix(fit, variable = "p_global")
  expect_snapshot(error = TRUE, extract_array_draws(dm, "p_global", dim1 = 5))
})
