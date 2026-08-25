mixture_data <- data.frame(
  d13C = c(-20, -21, -19, -22, -20.5, -21.5),
  d15N = c(10, 11, 9, 12, 10.5, 11.5),
  Region = factor(rep(c("A", "B"), each = 3))
)
source_data <- data.frame(
  Source = c("Beaver", "Deer"),
  d13C_mean = c(-25, -18), d13C_sd = c(1, 1),
  d15N_mean = c(5, 8), d15N_sd = c(1, 1)
)
tdf_data <- data.frame(
  Source = c("Beaver", "Deer"),
  d13C_mean = c(1, 1.2), d13C_sd = c(0.2, 0.3),
  d15N_mean = c(3, 3.1), d15N_sd = c(0.4, 0.5)
)

skip_if_not_installed("cmdstanr")
skip_if_not_installed("loo")

fit1 <- bsimm(
  formula = ~1, mixture_data = mixture_data, source_data = source_data,
  tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), source_means_sds = TRUE,
  chains = 1, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0,
  show_messages = FALSE, show_exceptions = FALSE
)
fit2 <- bsimm(
  formula = ~ 1 + (1 | Region), mixture_data = mixture_data, source_data = source_data,
  tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), source_means_sds = TRUE,
  chains = 1, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0,
  show_messages = FALSE, show_exceptions = FALSE
)

test_that("loo returns a psis_loo object", {
  l <- suppressWarnings(loo::loo(fit1))
  expect_true(loo::is.psis_loo(l))
})

test_that("waic returns a waic object", {
  w <- suppressWarnings(loo::waic(fit1))
  expect_true(loo::is.waic(w))
})

test_that("loo_compare ranks bsimms_fit objects by loo", {
  cmp <- suppressWarnings(loo::loo_compare(fit1, fit2, criterion = "loo"))
  expect_setequal(cmp$model, c("fit1", "fit2"))
})

test_that("loo_compare ranks bsimms_fit objects by waic", {
  cmp <- suppressWarnings(loo::loo_compare(fit1, fit2, criterion = "waic"))
  expect_setequal(cmp$model, c("fit1", "fit2"))
})

test_that("loo_compare accepts precomputed loo/waic objects", {
  l2 <- suppressWarnings(loo::loo(fit2))
  cmp <- suppressWarnings(
    loo::loo_compare(fit1, l2, criterion = "loo", model_names = c("fit1", "precomputed"))
  )
  expect_setequal(cmp$model, c("fit1", "precomputed"))
})

test_that("loo_compare errors when a precomputed object doesn't match criterion", {
  l1 <- suppressWarnings(loo::loo(fit1))
  expect_snapshot(error = TRUE, loo::loo_compare(fit2, l1, criterion = "waic"))
})

test_that("loo_compare respects custom model_names", {
  cmp <- suppressWarnings(loo::loo_compare(fit1, fit2, model_names = c("intercept_only", "region_re")))
  expect_setequal(cmp$model, c("intercept_only", "region_re"))
})

test_that("add_criterion caches loo by default", {
  fit1c <- suppressWarnings(add_criterion(fit1))
  expect_true(loo::is.psis_loo(fit1c$criteria$loo))
})

test_that("add_criterion can cache waic", {
  fit1c <- suppressWarnings(add_criterion(fit1, criterion = "waic"))
  expect_true(loo::is.waic(fit1c$criteria$waic))
})

test_that("add_criterion can cache both loo and waic at once", {
  fit1c <- suppressWarnings(add_criterion(fit1, criterion = c("loo", "waic")))
  expect_true(loo::is.psis_loo(fit1c$criteria$loo))
  expect_true(loo::is.waic(fit1c$criteria$waic))
})

test_that("add_criterion can cache bayes_R2 as raw draws", {
  skip_if_not_installed("rstantools")
  fit1c <- add_criterion(fit1, criterion = "bayes_R2")
  expect_true(is.matrix(fit1c$criteria$bayes_R2))
  expect_equal(colnames(fit1c$criteria$bayes_R2), c("d13C", "d15N"))
})

test_that("loo_compare informs when a criterion is computed on the fly", {
  expect_message(
    suppressWarnings(loo::loo_compare(fit1, fit2, criterion = "loo")),
    "not cached"
  )
})

test_that("loo_compare reuses a cached criterion without recomputing", {
  fit1c <- suppressWarnings(add_criterion(fit1, criterion = "loo"))
  msgs <- testthat::capture_messages(
    suppressWarnings(loo::loo_compare(fit1c, fit2, criterion = "loo"))
  )
  expect_true(any(grepl("fit2", msgs)))
  expect_false(any(grepl("fit1c", msgs)))
})
