mixture_data <- data.frame(
  d13C = c(-20, -21, -19, -22, -20.5, -21.5),
  d15N = c(10, 11, 9, 12, 10.5, 11.5),
  Region = factor(rep(c("A", "B"), each = 3)),
  Sex = factor(rep(c("M", "F"), 3))
)
source_data <- data.frame(
  Source = rep(c("Beaver", "Deer"), each = 3),
  d13C = c(-25, -24, -26, -18, -17, -19),
  d15N = c(5, 6, 4, 8, 9, 7),
  d13C_conc = c(0.3, 0.32, 0.28, 0.35, 0.33, 0.37)
)
tdf_data <- data.frame(
  Source = c("Beaver", "Deer"),
  d13C_mean = c(1, 1.2),
  d13C_sd = c(0.2, 0.3),
  d15N_mean = c(3, 3.1),
  d15N_sd = c(0.4, 0.5)
)

test_that("build_bsimms_spec derives dimensions and matrices correctly", {
  spec <- build_bsimms_spec(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data,
    tdf_data = tdf_data,
    isotope_names = c("d13C", "d15N")
  )
  expect_s3_class(spec, "bsimms_spec")
  expect_equal(spec$K, 2)
  expect_equal(spec$J, 2)
  expect_equal(spec$D, 1)
  expect_equal(spec$N, 6)
  expect_equal(spec$source_names, c("Beaver", "Deer"))
  expect_equal(dim(spec$V), c(2, 1))
  expect_equal(spec$source$mode, "raw")
  expect_equal(spec$tdf$mode, "summary")
})

test_that("the intercept is dropped from the fixed-effect design", {
  spec <- build_bsimms_spec(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data,
    tdf_data = tdf_data,
    isotope_names = c("d13C", "d15N")
  )
  expect_equal(spec$P, 0)
  expect_equal(ncol(spec$X), 0)
  expect_false("(Intercept)" %in% spec$fixed_names)
})

test_that("fixed covariates and group-level terms are captured", {
  spec <- build_bsimms_spec(
    formula = ~ Sex + (1 | Region),
    mixture_data = mixture_data,
    source_data = source_data,
    tdf_data = tdf_data,
    isotope_names = c("d13C", "d15N")
  )
  expect_true(any(grepl("^Sex", spec$fixed_names)))
  expect_equal(spec$P, ncol(spec$X))
  expect_length(spec$re_terms, 1)
  expect_equal(spec$re_terms[[1]]$group, "Region")
})

test_that("has_conc_dep reflects the conc_dep argument", {
  spec_off <- build_bsimms_spec(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data,
    tdf_data = tdf_data,
    isotope_names = "d13C",
    conc_dep = FALSE
  )
  spec_on <- build_bsimms_spec(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data,
    tdf_data = tdf_data,
    isotope_names = "d13C",
    conc_dep = TRUE
  )
  expect_false(spec_off$has_conc_dep)
  expect_null(spec_off$conc)
  expect_true(spec_on$has_conc_dep)
  expect_equal(dim(spec_on$conc), c(2L, 1L))
})

test_that("center/scale default to TRUE and are captured on the spec", {
  spec <- build_bsimms_spec(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data,
    tdf_data = tdf_data,
    isotope_names = c("d13C", "d15N")
  )
  expect_true(spec$center)
  expect_true(spec$scale)

  spec_off <- build_bsimms_spec(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data,
    tdf_data = tdf_data,
    isotope_names = c("d13C", "d15N"),
    center = FALSE,
    scale = FALSE
  )
  expect_false(spec_off$center)
  expect_false(spec_off$scale)

  expect_snapshot(
    error = TRUE,
    build_bsimms_spec(
      formula = ~1,
      mixture_data = mixture_data,
      source_data = source_data,
      tdf_data = tdf_data,
      isotope_names = c("d13C", "d15N"),
      center = "yes"
    )
  )
  expect_snapshot(
    error = TRUE,
    build_bsimms_spec(
      formula = ~1,
      mixture_data = mixture_data,
      source_data = source_data,
      tdf_data = tdf_data,
      isotope_names = c("d13C", "d15N"),
      scale = "yes"
    )
  )
})

test_that("a group-level term with a length mismatch against the mixture data is rejected", {
  expect_snapshot(
    error = TRUE,
    build_bsimms_spec(
      formula = ~ (1 | Region[1:3]),
      mixture_data = mixture_data,
      source_data = source_data,
      tdf_data = tdf_data,
      isotope_names = c("d13C", "d15N")
    )
  )
})

test_that("x_continuous is aligned with fixed_names and the dropped intercept", {
  mixture_data_x <- cbind(mixture_data, Age = c(1, 2, 3, 4, 5, 6))
  spec <- build_bsimms_spec(
    formula = ~ Age + Sex,
    mixture_data = mixture_data_x,
    source_data = source_data,
    tdf_data = tdf_data,
    isotope_names = c("d13C", "d15N")
  )
  expect_length(spec$x_continuous, spec$P)
  expect_true(spec$x_continuous[spec$fixed_names == "Age"])
  expect_false(any(spec$x_continuous[grepl("^Sex", spec$fixed_names)]))
})

test_that("has_beta_raw/has_ilr_global_raw reflect center/scale and x_continuous", {
  mixture_data_x <- cbind(mixture_data, Age = c(1, 2, 3, 4, 5, 6))
  spec_cont <- build_bsimms_spec(
    formula = ~Age,
    mixture_data = mixture_data_x,
    source_data = source_data,
    tdf_data = tdf_data,
    isotope_names = c("d13C", "d15N")
  )
  expect_true(spec_cont$has_beta_raw)
  expect_true(spec_cont$has_ilr_global_raw)

  spec_center_only <- build_bsimms_spec(
    formula = ~Age,
    mixture_data = mixture_data_x,
    source_data = source_data,
    tdf_data = tdf_data,
    isotope_names = c("d13C", "d15N"),
    scale = FALSE
  )
  expect_false(spec_center_only$has_beta_raw)
  expect_true(spec_center_only$has_ilr_global_raw)

  # requested but with nothing continuous to apply to: no effect
  spec_factor_only <- build_bsimms_spec(
    formula = ~Sex,
    mixture_data = mixture_data_x,
    source_data = source_data,
    tdf_data = tdf_data,
    isotope_names = c("d13C", "d15N")
  )
  expect_false(spec_factor_only$has_beta_raw)
  expect_false(spec_factor_only$has_ilr_global_raw)
})

test_that("error_structure defaults to process_residual and is validated", {
  spec <- build_bsimms_spec(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data,
    tdf_data = tdf_data,
    isotope_names = c("d13C", "d15N")
  )
  expect_equal(spec$error_structure, "process_residual")
  expect_snapshot(
    error = TRUE,
    build_bsimms_spec(
      formula = ~1,
      mixture_data = mixture_data,
      source_data = source_data,
      tdf_data = tdf_data,
      isotope_names = c("d13C", "d15N"),
      error_structure = "banana"
    )
  )
})

test_that("a single mixture data point requires process_only error", {
  one_mix <- mixture_data[1, ]
  spec <- build_bsimms_spec(
    formula = ~1,
    mixture_data = one_mix,
    source_data = source_data,
    tdf_data = tdf_data,
    isotope_names = c("d13C", "d15N"),
    error_structure = "process_only"
  )
  expect_equal(spec$N, 1)
  expect_snapshot(
    error = TRUE,
    build_bsimms_spec(
      formula = ~1,
      mixture_data = one_mix,
      source_data = source_data,
      tdf_data = tdf_data,
      isotope_names = c("d13C", "d15N"),
      error_structure = "process_residual"
    )
  )
  expect_snapshot(
    error = TRUE,
    build_bsimms_spec(
      formula = ~1,
      mixture_data = one_mix,
      source_data = source_data,
      tdf_data = tdf_data,
      isotope_names = c("d13C", "d15N"),
      error_structure = "residual_only"
    )
  )
})

test_that("one mixture sample per level of a single fixed factor requires process_only", {
  ind_mix <- data.frame(
    d13C = c(-20, -21, -19),
    d15N = c(10, 11, 9),
    Individual = factor(c("A", "B", "C"))
  )
  spec <- build_bsimms_spec(
    formula = ~Individual,
    mixture_data = ind_mix,
    source_data = source_data,
    tdf_data = tdf_data,
    isotope_names = c("d13C", "d15N"),
    error_structure = "process_only"
  )
  expect_equal(spec$N, 3)
  expect_snapshot(
    error = TRUE,
    build_bsimms_spec(
      formula = ~Individual,
      mixture_data = ind_mix,
      source_data = source_data,
      tdf_data = tdf_data,
      isotope_names = c("d13C", "d15N"),
      error_structure = "process_residual"
    )
  )
})

test_that("one mixture sample per level of a random-effect group requires process_only", {
  ind_mix <- mixture_data
  ind_mix$Individual <- factor(paste0("I", seq_len(nrow(ind_mix))))
  spec <- build_bsimms_spec(
    formula = ~ 1 + (1 | Individual),
    mixture_data = ind_mix,
    source_data = source_data,
    tdf_data = tdf_data,
    isotope_names = c("d13C", "d15N"),
    error_structure = "process_only"
  )
  expect_equal(spec$N, 6)
  expect_snapshot(
    error = TRUE,
    build_bsimms_spec(
      formula = ~ 1 + (1 | Individual),
      mixture_data = ind_mix,
      source_data = source_data,
      tdf_data = tdf_data,
      isotope_names = c("d13C", "d15N"),
      error_structure = "process_residual"
    )
  )
})

test_that("build_bsimms_spec errors on invalid isotope_names", {
  expect_snapshot(
    error = TRUE,
    build_bsimms_spec(
      formula = ~1,
      mixture_data = mixture_data,
      source_data = source_data,
      tdf_data = tdf_data,
      isotope_names = character(0)
    )
  )
})

test_that("build_bsimms_spec errors with fewer than 2 sources", {
  one_source <- source_data[source_data$Source == "Beaver", ]
  one_tdf <- tdf_data[tdf_data$Source == "Beaver", ]
  expect_snapshot(
    error = TRUE,
    build_bsimms_spec(
      formula = ~1,
      mixture_data = mixture_data,
      source_data = one_source,
      tdf_data = one_tdf,
      isotope_names = c("d13C", "d15N")
    )
  )
})

test_that("build_bsimms_spec errors when the formula's design matrix drops rows", {
  d <- mixture_data
  d$Region[1] <- NA
  expect_snapshot(
    error = TRUE,
    build_bsimms_spec(
      formula = ~Region,
      mixture_data = d,
      source_data = source_data,
      tdf_data = tdf_data,
      isotope_names = c("d13C", "d15N")
    )
  )
})
