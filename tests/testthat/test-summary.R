mixture_data <- data.frame(
  d13C = c(-20, -21, -19, -22, -20.5, -21.5),
  d15N = c(10, 11, 9, 12, 10.5, 11.5),
  elevation = c(100, 120, 110, 300, 320, 310),
  Region = factor(rep(c("A", "B"), each = 3))
)
source_data <- data.frame(
  Source = rep(c("Beaver", "Deer", "Otter"), each = 3),
  d13C = c(-25, -24, -26, -18, -17, -19, -21, -20, -22),
  d15N = c(5, 6, 4, 8, 9, 7, 6, 7, 5)
)
tdf_data <- data.frame(
  Source = c("Beaver", "Deer", "Otter"),
  d13C_mean = c(1, 1.2, 1.1), d13C_sd = c(0.2, 0.3, 0.25),
  d15N_mean = c(3, 3.1, 3.05), d15N_sd = c(0.4, 0.5, 0.45)
)

skip_if_not_installed("cmdstanr")

fit_full <- bsimm(
  formula = ~ elevation + (1 | Region), mixture_data = mixture_data, source_data = source_data,
  tdf_data = tdf_data, isotope_names = c("d13C", "d15N"),
  chains = 1, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0
)
fit_min <- bsimm(
  formula = ~1, mixture_data = mixture_data, source_data = source_data,
  tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), error_structure = "process_only",
  chains = 1, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0
)

test_that("print.bsimms_fit prints a one-screen model overview", {
  expect_snapshot(print(fit_full))
})

test_that("summary.bsimms_fit includes fixed/group/error sections when present", {
  s <- summary(fit_full)
  expect_s3_class(s, "summary.bsimms_fit")
  expect_equal(s$population_proportions$source, c("Beaver", "Deer", "Otter"))
  expect_equal(s$fixed$fixed_effect, rep("elevation", 2))
  expect_equal(s$fixed$ilr_dim, c(1, 2))
  expect_equal(s$group$group, rep("Region", 2))
  expect_equal(s$group$term, rep("(Intercept)", 2))
  expect_equal(s$group$ilr_dim, c(1, 2))
  expect_null(s$group$variable)
  expect_equal(s$error$isotope, c("d13C", "d15N"))
})

test_that("summary.bsimms_fit omits fixed/group/error sections when absent", {
  s <- summary(fit_min)
  expect_null(s$fixed)
  expect_null(s$group)
  expect_null(s$error)
  expect_equal(s$population_proportions$source, c("Beaver", "Deer", "Otter"))
})

test_that("summary.bsimms_fit's robust argument switches mean/sd for median/mad", {
  s_mean <- summary(fit_full, robust = FALSE)
  s_median <- summary(fit_full, robust = TRUE)
  expect_true(all(c("mean", "sd") %in% names(s_mean$population_proportions)))
  expect_true(all(c("median", "mad") %in% names(s_median$population_proportions)))
})

test_that("summary.bsimms_fit's probs argument controls the quantile columns", {
  s <- summary(fit_full, probs = c(0.1, 0.9))
  expect_true(all(c("q10", "q90") %in% names(s$population_proportions)))
})

test_that("print.summary.bsimms_fit prints every section present in the summary", {
  out <- capture.output(print(summary(fit_full)))
  expect_true(any(grepl("^Population-average source proportions:", out)))
  expect_true(any(grepl("^Fixed effects", out)))
  expect_true(any(grepl("^Group-level standard deviations", out)))
  expect_true(any(grepl("^Error term\\(s\\):", out)))
  group_header <- out[which(grepl("^Group-level standard deviations", out)) + 1]
  expect_true(grepl("group", group_header))
  expect_true(grepl("term", group_header))
  expect_true(grepl("ilr_dim", group_header))
})

test_that("print.summary.bsimms_fit omits sections absent from the summary", {
  out <- capture.output(print(summary(fit_min)))
  expect_true(any(grepl("^Population-average source proportions:", out)))
  expect_false(any(grepl("^Fixed effects", out)))
  expect_false(any(grepl("^Group-level standard deviations", out)))
  expect_false(any(grepl("^Error term\\(s\\):", out)))
})

test_that("print.summary.bsimms_fit rounds to the requested number of digits", {
  out <- capture.output(print(summary(fit_full), digits = 1))
  pp_line <- out[grepl("^ Beaver|^ Deer|^ Otter", out)]
  nums <- regmatches(pp_line, gregexpr("-?[0-9]+\\.[0-9]+", pp_line))
  decimals <- vapply(unlist(nums), function(n) nchar(strsplit(n, "\\.")[[1]][2]), integer(1))
  expect_true(all(decimals <= 1))
})
