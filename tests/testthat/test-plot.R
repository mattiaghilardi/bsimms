mixture_data <- data.frame(
  d13C = c(-20, -21, -19, -22, -20.5, -21.5),
  d15N = c(10, 11, 9, 12, 10.5, 11.5),
  Region = factor(rep(c("A", "B"), each = 3))
)
source_data <- data.frame(
  Source = c("Beaver", "Deer"),
  d13C_mean = c(-25, -18),
  d13C_sd = c(1, 1),
  d15N_mean = c(5, 8),
  d15N_sd = c(1, 1)
)
tdf_data <- data.frame(
  Source = c("Beaver", "Deer"),
  d13C_mean = c(1, 1.2),
  d13C_sd = c(0.2, 0.3),
  d15N_mean = c(3, 3.1),
  d15N_sd = c(0.4, 0.5)
)

skip_if_not_installed("cmdstanr")
skip_if_not_installed("bayesplot")

grDevices::pdf(nullfile())

fit <- bsimm(
  formula = ~ 1 + (1 | Region),
  mixture_data = mixture_data,
  source_data = source_data,
  tdf_data = tdf_data,
  isotope_names = c("d13C", "d15N"),
  source_means_sds = TRUE,
  chains = 1,
  iter_warmup = 200,
  iter_sampling = 100,
  seed = 1,
  refresh = 0,
  show_messages = FALSE,
  show_exceptions = FALSE
)
fit_1iso <- bsimm(
  formula = ~1,
  mixture_data = mixture_data,
  source_data = source_data,
  tdf_data = tdf_data,
  isotope_names = "d13C",
  source_means_sds = TRUE,
  chains = 1,
  iter_warmup = 200,
  iter_sampling = 100,
  seed = 1,
  refresh = 0,
  show_messages = FALSE,
  show_exceptions = FALSE
)

# plot.bsimms_fit ---------------------------------------------------------

test_that("plot.bsimms_fit builds one page for the default variable set", {
  plots <- plot(fit, plot = FALSE)
  expect_type(plots, "list")
  expect_length(plots, 1)
  expect_s3_class(plots[[1]], "bayesplot_grid")
})

test_that("plot.bsimms_fit respects an explicit variable selection", {
  plots <- plot(fit, variable = "p_global", plot = FALSE)
  expect_length(plots, 1)
})

test_that("plot.bsimms_fit paginates when nvariables is smaller than the parameter count", {
  plots <- plot(fit, variable = "p_global", nvariables = 1, plot = FALSE)
  expect_length(plots, 2)
})

test_that("plot.bsimms_fit errors when no parameters match", {
  expect_snapshot(error = TRUE, plot(fit, variable = "banana", plot = FALSE))
})

test_that("plot.bsimms_fit can actually draw to a graphics device", {
  grDevices::pdf(nullfile())
  on.exit(grDevices::dev.off())
  plots <- plot(fit, variable = "p_global", plot = TRUE, ask = FALSE)
  expect_length(plots, 1)
})

# pp_check.bsimms_fit -------------------------------------------------------

test_that("pp_check.bsimms_fit requires resp for multi-isotope models", {
  expect_snapshot(error = TRUE, bayesplot::pp_check(fit))
})

test_that("pp_check.bsimms_fit defaults resp for single-isotope models", {
  g <- bayesplot::pp_check(fit_1iso)
  expect_s3_class(g, "ggplot")
})

test_that("pp_check.bsimms_fit works with an explicit resp and type", {
  g <- bayesplot::pp_check(fit, resp = "d15N", type = "hist")
  expect_s3_class(g, "ggplot")
})

test_that("pp_check.bsimms_fit errors on an invalid type", {
  expect_snapshot(
    error = TRUE,
    bayesplot::pp_check(fit, resp = "d13C", type = "banana")
  )
})

test_that("pp_check.bsimms_fit subsamples with ndraws", {
  g <- bayesplot::pp_check(fit_1iso, ndraws = 10)
  expect_s3_class(g, "ggplot")
})

test_that("pp_check.bsimms_fit defaults to 10 draws for overlay-style types", {
  expect_message(
    bayesplot::pp_check(fit_1iso, type = "hist"),
    "Using 10 posterior draws"
  )
})

test_that("pp_check.bsimms_fit defaults to all draws for aggregate types", {
  expect_message(
    bayesplot::pp_check(fit_1iso, type = "stat"),
    "Using all posterior draws"
  )
})

test_that("pp_check.bsimms_fit doesn't message when ndraws is supplied explicitly", {
  expect_no_message(bayesplot::pp_check(fit_1iso, ndraws = 5))
})

test_that("pp_check.bsimms_fit errors when ndraws exceeds the number of draws available", {
  expect_snapshot(error = TRUE, bayesplot::pp_check(fit_1iso, ndraws = 1000))
})

# plot_proportions ------------------------------------------------------------

skip_if_not_installed("ggplot2")

p_arr <- posterior_proportions(fit)

test_that("plot_proportions errors when p_arr lacks source-name dimnames", {
  bare_arr <- p_arr
  dimnames(bare_arr) <- NULL
  expect_snapshot(error = TRUE, plot_proportions(bare_arr))
})

test_that("plot_proportions builds a density plot for a single observation", {
  g <- plot_proportions(p_arr[, 1, , drop = FALSE], type = "density")
  expect_s3_class(g, "ggplot")
})

test_that("plot_proportions builds a histogram for a single observation", {
  g <- plot_proportions(p_arr[, 1, , drop = FALSE], type = "histogram")
  expect_s3_class(g, "ggplot")
})

test_that("plot_proportions errors for density/histogram with more than one observation", {
  expect_snapshot(error = TRUE, plot_proportions(p_arr, type = "density"))
})

test_that("plot_proportions builds an interval plot for multiple observations", {
  g <- plot_proportions(p_arr, type = "interval")
  expect_s3_class(g, "ggplot")
})

test_that("plot_proportions builds an interval plot for a single observation too", {
  g <- plot_proportions(p_arr[, 1, , drop = FALSE], type = "interval")
  expect_s3_class(g, "ggplot")
})

test_that("plot_proportions errors on invalid probs", {
  expect_snapshot(
    error = TRUE,
    plot_proportions(p_arr, type = "interval", probs = c(0, 0.95))
  )
})

test_that("plot_proportions's ... is forwarded to geom_histogram", {
  g <- plot_proportions(
    p_arr[, 1, , drop = FALSE],
    type = "histogram",
    bins = 10
  )
  expect_equal(g$layers[[1]]$stat_params$bins, 10)
})

test_that("plot_proportions's ... is forwarded to geom_linerange", {
  g <- plot_proportions(p_arr, type = "interval", linetype = "dashed")
  expect_equal(g$layers[[1]]$aes_params$linetype, "dashed")
})
