test_that("default call returns well-formed data and truth", {
  sim <- simulate_bsimms_data(n_mixture_obs = 20, seed = 1)
  expect_named(sim, c(
    "formula", "mixture_data", "source_data", "tdf_data", "isotope_names",
    "source_means_sds", "tdf_means_sds", "conc_dep", "error_structure", "source_col", "truth"
  ))
  expect_equal(nrow(sim$mixture_data), 20)
  expect_equal(sim$isotope_names, c("isotope1", "isotope2"))
  expect_setequal(unique(sim$source_data$Source), c("source1", "source2", "source3"))
  expect_null(sim$truth$fixed)
  expect_null(sim$truth$random)
  expect_length(sim$truth$p_global, 3)
  expect_equal(sum(sim$truth$p_global), 1)
})

test_that("source_names/isotope_names override n_sources/n_isotopes", {
  sim <- simulate_bsimms_data(
    n_mixture_obs = 10, source_names = c("Beaver", "Deer"),
    isotope_names = c("d13C", "d15N"), seed = 1
  )
  expect_equal(sim$isotope_names, c("d13C", "d15N"))
  expect_setequal(unique(sim$source_data$Source), c("Beaver", "Deer"))
})

test_that("fixed and random-effect covariates are generated correctly", {
  sim <- simulate_bsimms_data(
    ~ elevation + Sex + (1 | Region),
    n_mixture_obs = 60, n_levels = list(Sex = 2), n_groups = list(Region = 3), seed = 1
  )
  expect_true(is.numeric(sim$mixture_data$elevation))
  expect_s3_class(sim$mixture_data$Sex, "factor")
  expect_equal(nlevels(sim$mixture_data$Sex), 2)
  expect_s3_class(sim$mixture_data$Region, "factor")
  expect_equal(nlevels(sim$mixture_data$Region), 3)
  expect_equal(dim(sim$truth$fixed), c(2, 2))
  expect_length(sim$truth$random, 1)
  expect_equal(sim$truth$random[[1]]$group, "Region")
  expect_equal(dim(sim$truth$random[[1]]$effects), c(2, 3))
})

test_that("balanced levels are split as evenly as possible", {
  sim <- simulate_bsimms_data(~Region, n_mixture_obs = 11, n_levels = list(Region = 3), seed = 1)
  expect_equal(sort(as.integer(table(sim$mixture_data$Region))), c(3, 4, 4))
})

test_that("unbalanced levels still get at least one observation each", {
  sim <- simulate_bsimms_data(
    ~Region,
    n_mixture_obs = 6, n_levels = list(Region = 5), balanced = FALSE, seed = 1
  )
  counts <- table(sim$mixture_data$Region)
  expect_equal(length(counts), 5)
  expect_true(all(counts >= 1))
  expect_equal(sum(counts), 6)
})

test_that("unknown n_levels/n_groups names are rejected", {
  expect_snapshot(error = TRUE, simulate_bsimms_data(~Sex, n_mixture_obs = 10, n_levels = list(banana = 2)))
  expect_snapshot(
    error = TRUE,
    simulate_bsimms_data(~ (1 | Region), n_mixture_obs = 10, n_groups = list(banana = 2))
  )
})

test_that("a grouping factor missing from n_groups is rejected", {
  expect_snapshot(error = TRUE, simulate_bsimms_data(~ (1 | Region), n_mixture_obs = 10))
})

test_that("random slopes are rejected", {
  expect_snapshot(
    error = TRUE,
    simulate_bsimms_data(~ (elevation | Region), n_mixture_obs = 10, n_groups = list(Region = 2))
  )
})

test_that("n_levels/n_groups entries must be integers >= 2", {
  expect_snapshot(error = TRUE, simulate_bsimms_data(~Sex, n_mixture_obs = 10, n_levels = list(Sex = 1)))
})

test_that("too many levels for the available observations is rejected", {
  expect_snapshot(error = TRUE, simulate_bsimms_data(~Region, n_mixture_obs = 2, n_levels = list(Region = 5)))
})

test_that("p_global can be overridden and is validated", {
  sim <- simulate_bsimms_data(n_mixture_obs = 10, n_sources = 3, p_global = c(0.2, 0.5, 0.3), seed = 1)
  expect_equal(sim$truth$p_global, c(0.2, 0.5, 0.3))
  expect_snapshot(
    error = TRUE,
    simulate_bsimms_data(n_mixture_obs = 10, n_sources = 3, p_global = c(0.2, 0.5, 0.5))
  )
})

test_that("sigma/resid_prop truth is only present for the relevant error_structure", {
  sim_resid <- simulate_bsimms_data(n_mixture_obs = 10, error_structure = "residual_only", seed = 1)
  expect_named(sim_resid$truth$sigma, c("isotope1", "isotope2"))
  expect_null(sim_resid$truth$resid_prop)

  sim_proc <- simulate_bsimms_data(n_mixture_obs = 10, error_structure = "process_only", seed = 1)
  expect_null(sim_proc$truth$sigma)
  expect_null(sim_proc$truth$resid_prop)

  sim_pr <- simulate_bsimms_data(n_mixture_obs = 10, error_structure = "process_residual", seed = 1)
  expect_named(sim_pr$truth$resid_prop, c("isotope1", "isotope2"))
  expect_null(sim_pr$truth$sigma)
})

test_that("source/tdf mean-sd truth is only kept when the corresponding data is raw", {
  sim <- simulate_bsimms_data(n_mixture_obs = 10, source_means_sds = TRUE, tdf_means_sds = FALSE, seed = 1)
  expect_null(sim$truth$source_mean)
  expect_null(sim$truth$source_sd)
  expect_equal(dim(sim$truth$tdf_mean), c(3, 2))
  expect_equal(dim(sim$truth$tdf_sd), c(3, 2))
})

test_that("raw vs summarised source/tdf data have the expected shape", {
  sim_raw <- simulate_bsimms_data(n_mixture_obs = 10, n_sources = 2, n_source_obs = 7, source_means_sds = FALSE, seed = 1)
  expect_equal(nrow(sim_raw$source_data), 2 * 7)
  expect_named(sim_raw$source_data, c("Source", "isotope1", "isotope2"))

  sim_summary <- simulate_bsimms_data(n_mixture_obs = 10, n_sources = 2, source_means_sds = TRUE, seed = 1)
  expect_equal(nrow(sim_summary$source_data), 2)
  expect_named(sim_summary$source_data, c("Source", "isotope1_mean", "isotope1_sd", "isotope2_mean", "isotope2_sd"))
})

test_that("n_source_obs/n_tdf_obs allow unbalanced per-source replicate counts", {
  sim <- simulate_bsimms_data(
    n_mixture_obs = 10, n_sources = 3, n_source_obs = c(5, 20, 8),
    tdf_means_sds = FALSE, n_tdf_obs = c(3, 4, 5), seed = 1
  )
  expect_equal(as.integer(table(sim$source_data$Source)), c(5, 20, 8))
  expect_equal(as.integer(table(sim$tdf_data$Source)), c(3, 4, 5))
})

test_that("n_source_obs/n_tdf_obs reject the wrong length or non-positive values", {
  expect_snapshot(
    error = TRUE,
    simulate_bsimms_data(n_mixture_obs = 10, n_sources = 3, n_source_obs = c(5, 20))
  )
  expect_snapshot(
    error = TRUE,
    simulate_bsimms_data(n_mixture_obs = 10, n_sources = 3, n_source_obs = c(5, -1, 8))
  )
})

test_that("conc_dep adds <isotope>_conc columns for raw and summarised source data", {
  sim_raw <- simulate_bsimms_data(n_mixture_obs = 10, conc_dep = TRUE, source_means_sds = FALSE, seed = 1)
  expect_true(all(c("isotope1_conc", "isotope2_conc") %in% names(sim_raw$source_data)))
  expect_true(all(sim_raw$source_data$isotope1_conc > 0 & sim_raw$source_data$isotope1_conc <= 1))

  sim_summary <- simulate_bsimms_data(n_mixture_obs = 10, conc_dep = TRUE, source_means_sds = TRUE, seed = 1)
  expect_true(all(c("isotope1_conc", "isotope2_conc") %in% names(sim_summary$source_data)))
})

test_that("an underdetermined system (n_sources > n_isotopes + 1) warns", {
  expect_snapshot(invisible(simulate_bsimms_data(n_mixture_obs = 10, n_sources = 4, n_isotopes = 2, seed = 1)))
  expect_no_warning(simulate_bsimms_data(n_mixture_obs = 10, n_sources = 3, n_isotopes = 2, seed = 1))
})

test_that("source means are separated and non-collinear for K >= 3, J >= 2", {
  sim <- simulate_bsimms_data(n_mixture_obs = 10, n_sources = 4, n_isotopes = 3, seed = 1)
  m <- sim$truth$source_mean
  expect_gt(min(svd(scale(m, scale = FALSE))$d), 0)
})

test_that("source_col is honoured and echoed back", {
  sim <- simulate_bsimms_data(n_mixture_obs = 10, n_sources = 2, source_col = "SourceID", seed = 1)
  expect_equal(sim$source_col, "SourceID")
  expect_true("SourceID" %in% names(sim$source_data))
  expect_true("SourceID" %in% names(sim$tdf_data))
})

test_that("seed makes output reproducible and does not leak into the global RNG state", {
  sim1 <- simulate_bsimms_data(n_mixture_obs = 10, seed = 42)
  sim2 <- simulate_bsimms_data(n_mixture_obs = 10, seed = 42)
  expect_equal(sim1$mixture_data, sim2$mixture_data)

  set.seed(1)
  before <- .Random.seed
  simulate_bsimms_data(n_mixture_obs = 10, seed = 99)
  expect_identical(.Random.seed, before)
})

test_that("simulated data validates against make_standata() across configurations", {
  configs <- expand.grid(
    source_means_sds = c(FALSE, TRUE),
    tdf_means_sds = c(FALSE, TRUE),
    error_structure = c("residual_only", "process_only", "process_residual"),
    conc_dep = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(configs))) {
    cfg <- configs[i, ]
    sim <- simulate_bsimms_data(
      ~ elevation + (1 | Region),
      n_mixture_obs = 30, n_groups = list(Region = 3),
      source_means_sds = cfg$source_means_sds, tdf_means_sds = cfg$tdf_means_sds,
      error_structure = cfg$error_structure, conc_dep = cfg$conc_dep, seed = 1
    )
    sdata <- make_standata(
      sim$formula,
      mixture_data = sim$mixture_data, source_data = sim$source_data, tdf_data = sim$tdf_data,
      isotope_names = sim$isotope_names, source_means_sds = sim$source_means_sds,
      tdf_means_sds = sim$tdf_means_sds, error_structure = sim$error_structure, conc_dep = sim$conc_dep
    )
    expect_s3_class(sdata, "bsimms_standata")
  }
})

test_that("simulated data fits via bsimm()", {
  skip_if_not_installed("cmdstanr")
  sim <- simulate_bsimms_data(
    ~ elevation + Sex + (1 | Region),
    n_mixture_obs = 40, n_levels = list(Sex = 2), n_groups = list(Region = 3), seed = 1
  )
  fit <- bsimm(
    sim$formula,
    mixture_data = sim$mixture_data, source_data = sim$source_data, tdf_data = sim$tdf_data,
    isotope_names = sim$isotope_names, source_means_sds = sim$source_means_sds,
    tdf_means_sds = sim$tdf_means_sds, error_structure = sim$error_structure,
    chains = 1, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0
  )
  expect_s3_class(fit, "bsimms_fit")
})
