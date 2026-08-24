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
skip_if_not_installed("rstantools")

fit <- bsimm(
  formula = ~ 1 + (1 | Region), mixture_data = mixture_data, source_data = source_data,
  tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), source_means_sds = TRUE,
  chains = 1, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0,
  show_messages = FALSE, show_exceptions = FALSE
)

test_that("bayes_R2 returns a summary data frame for every isotope by default", {
  r2 <- rstantools::bayes_R2(fit)
  expect_s3_class(r2, "data.frame")
  expect_equal(r2$isotope, c("d13C", "d15N"))
  expect_true(all(c("mean", "sd", "q2.5", "q97.5") %in% names(r2)))
})

test_that("bayes_R2 respects an explicit resp", {
  r2 <- rstantools::bayes_R2(fit, resp = "d15N")
  expect_equal(r2$isotope, "d15N")
})

test_that("bayes_R2 errors on an invalid resp", {
  expect_snapshot(error = TRUE, rstantools::bayes_R2(fit, resp = "banana"))
})

test_that("bayes_R2's robust argument switches mean/sd for median/mad", {
  r2 <- rstantools::bayes_R2(fit, robust = TRUE)
  expect_true(all(c("median", "mad", "q2.5", "q97.5") %in% names(r2)))
})

test_that("bayes_R2's probs argument controls the quantile columns", {
  r2 <- rstantools::bayes_R2(fit, probs = c(0.1, 0.9))
  expect_true(all(c("q10", "q90") %in% names(r2)))
})

test_that("bayes_R2 returns raw draws when summary = FALSE", {
  r2 <- rstantools::bayes_R2(fit, summary = FALSE)
  expect_true(is.matrix(r2))
  expect_equal(colnames(r2), c("d13C", "d15N"))
  expect_equal(nrow(r2), 100)
  expect_true(all(r2 >= 0 & r2 <= 1))
})
