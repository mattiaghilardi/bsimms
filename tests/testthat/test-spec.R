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
  d13C_mean = c(1, 1.2), d13C_sd = c(0.2, 0.3),
  d15N_mean = c(3, 3.1), d15N_sd = c(0.4, 0.5)
)

test_that("build_bsimms_spec derives dimensions and matrices correctly", {
  spec <- build_bsimms_spec(
    formula = ~1, mixture_data = mixture_data, source_data = source_data,
    tdf_data = tdf_data, isotope_names = c("d13C", "d15N")
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
    formula = ~1, mixture_data = mixture_data, source_data = source_data,
    tdf_data = tdf_data, isotope_names = c("d13C", "d15N")
  )
  expect_equal(spec$P, 0)
  expect_equal(ncol(spec$X), 0)
  expect_false("(Intercept)" %in% spec$fixed_names)
})

test_that("fixed covariates and group-level terms are captured", {
  spec <- build_bsimms_spec(
    formula = ~ Sex + (1 | Region), mixture_data = mixture_data, source_data = source_data,
    tdf_data = tdf_data, isotope_names = c("d13C", "d15N")
  )
  expect_true(any(grepl("^Sex", spec$fixed_names)))
  expect_equal(spec$P, ncol(spec$X))
  expect_length(spec$re_terms, 1)
  expect_equal(spec$re_terms[[1]]$group, "Region")
})

test_that("has_conc_dep reflects the conc_dep argument", {
  spec_off <- build_bsimms_spec(
    formula = ~1, mixture_data = mixture_data, source_data = source_data,
    tdf_data = tdf_data, isotope_names = "d13C", conc_dep = FALSE
  )
  spec_on <- build_bsimms_spec(
    formula = ~1, mixture_data = mixture_data, source_data = source_data,
    tdf_data = tdf_data, isotope_names = "d13C", conc_dep = TRUE
  )
  expect_false(spec_off$has_conc_dep)
  expect_null(spec_off$conc)
  expect_true(spec_on$has_conc_dep)
  expect_equal(dim(spec_on$conc), c(2L, 1L))
})

test_that("error_structure defaults to process_residual and is validated", {
  spec <- build_bsimms_spec(
    formula = ~1, mixture_data = mixture_data, source_data = source_data,
    tdf_data = tdf_data, isotope_names = c("d13C", "d15N")
  )
  expect_equal(spec$error_structure, "process_residual")
  expect_snapshot(
    error = TRUE,
    build_bsimms_spec(
      formula = ~1, mixture_data = mixture_data, source_data = source_data,
      tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), error_structure = "banana"
    )
  )
})

test_that("build_bsimms_spec errors on invalid isotope_names", {
  expect_snapshot(
    error = TRUE,
    build_bsimms_spec(
      formula = ~1, mixture_data = mixture_data, source_data = source_data,
      tdf_data = tdf_data, isotope_names = character(0)
    )
  )
})

test_that("build_bsimms_spec errors with fewer than 2 sources", {
  one_source <- source_data[source_data$Source == "Beaver", ]
  one_tdf <- tdf_data[tdf_data$Source == "Beaver", ]
  expect_snapshot(
    error = TRUE,
    build_bsimms_spec(
      formula = ~1, mixture_data = mixture_data, source_data = one_source,
      tdf_data = one_tdf, isotope_names = c("d13C", "d15N")
    )
  )
})

test_that("build_bsimms_spec errors when the formula's design matrix drops rows", {
  d <- mixture_data
  d$Region[1] <- NA
  expect_snapshot(
    error = TRUE,
    build_bsimms_spec(
      formula = ~Region, mixture_data = d, source_data = source_data,
      tdf_data = tdf_data, isotope_names = c("d13C", "d15N")
    )
  )
})
