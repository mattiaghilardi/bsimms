# infer_source_names() ------------------------------------------------------

test_that("infer_source_names returns sorted unique source names", {
  sources <- data.frame(Source = c("Deer", "Beaver", "Deer", "Hare"))
  expect_identical(
    infer_source_names(source_data = sources, tdf_data = NULL),
    c("Beaver", "Deer", "Hare")
  )
})

test_that("infer_source_names cross-checks against tdf_data", {
  sources <- data.frame(Source = c("Beaver", "Deer"))
  tdf <- data.frame(Source = c("Deer", "Beaver"))
  expect_identical(
    infer_source_names(source_data = sources, tdf_data = tdf),
    c("Beaver", "Deer")
  )
})

test_that("infer_source_names errors on a non-string source_col", {
  sources <- data.frame(Source = c("Beaver", "Deer"))
  expect_snapshot(
    error = TRUE,
    infer_source_names(
      source_data = sources,
      tdf_data = NULL,
      source_col = c("Source", "Other")
    )
  )
  expect_snapshot(
    error = TRUE,
    infer_source_names(
      source_data = sources,
      tdf_data = NULL,
      source_col = NULL
    )
  )
})

test_that("infer_source_names errors on missing source_col", {
  expect_snapshot(
    error = TRUE,
    infer_source_names(source_data = data.frame(x = 1), tdf_data = NULL)
  )
  expect_snapshot(
    error = TRUE,
    infer_source_names(
      source_data = data.frame(Source = "A"),
      tdf_data = data.frame(x = 1)
    )
  )
})

test_that("infer_source_names errors on source/tdf name mismatch", {
  sources <- data.frame(Source = c("Beaver", "Deer"))
  tdf <- data.frame(Source = c("Beaver", "Hare"))
  expect_snapshot(
    error = TRUE,
    infer_source_names(source_data = sources, tdf_data = tdf)
  )
})

# prep_mixture_isotopes() ----------------------------------------------------

test_that("prep_mixture_isotopes builds a numeric matrix in isotope_names order", {
  d <- data.frame(d15N = c(1, 2), extra = c(9, 9), d13C = c(3, 4))
  y <- prep_mixture_isotopes(
    mixture_data = d,
    isotope_names = c("d13C", "d15N")
  )
  expect_equal(
    y,
    matrix(c(3, 4, 1, 2), nrow = 2, dimnames = list(NULL, c("d13C", "d15N")))
  )
})

test_that("prep_mixture_isotopes errors on missing isotope columns", {
  expect_snapshot(
    error = TRUE,
    prep_mixture_isotopes(
      mixture_data = data.frame(d13C = 1),
      isotope_names = c("d13C", "d15N")
    )
  )
})

test_that("prep_mixture_isotopes errors on missing values", {
  d <- data.frame(d13C = c(1, NA), d15N = c(1, 2))
  expect_snapshot(
    error = TRUE,
    prep_mixture_isotopes(mixture_data = d, isotope_names = c("d13C", "d15N"))
  )
})

# prep_iso_table(): raw mode -------------------------------------------------

test_that("prep_iso_table (raw) builds source_idx and Y correctly", {
  d <- data.frame(
    Source = c("Deer", "Beaver", "Deer"),
    d13C = c(1, 2, 3),
    d15N = c(4, 5, 6)
  )
  out <- prep_iso_table(
    data = d,
    isotope_names = c("d13C", "d15N"),
    source_names = c("Beaver", "Deer"),
    means_sds = FALSE
  )
  expect_equal(out$mode, "raw")
  expect_equal(out$n, 3)
  expect_equal(out$source_idx, c(2, 1, 2))
  expect_equal(unname(out$Y), matrix(c(1, 2, 3, 4, 5, 6), nrow = 3))
})

test_that("prep_iso_table errors on an invalid label", {
  d <- data.frame(Source = "Beaver", d13C = 1)
  expect_snapshot(
    error = TRUE,
    prep_iso_table(
      data = d,
      isotope_names = "d13C",
      source_names = "Beaver",
      means_sds = FALSE,
      label = "banana"
    )
  )
})

test_that("prep_iso_table (raw) errors on missing source_col", {
  expect_snapshot(
    error = TRUE,
    prep_iso_table(
      data = data.frame(x = 1),
      isotope_names = "d13C",
      source_names = "Beaver",
      means_sds = FALSE
    )
  )
})

test_that("prep_iso_table (raw) errors on unknown source", {
  d <- data.frame(Source = "Fox", d13C = 1)
  expect_snapshot(
    error = TRUE,
    prep_iso_table(
      data = d,
      isotope_names = "d13C",
      source_names = c("Beaver", "Deer"),
      means_sds = FALSE
    )
  )
})

test_that("prep_iso_table (raw) errors on missing isotope columns", {
  d <- data.frame(Source = "Beaver", d13C = 1)
  expect_snapshot(
    error = TRUE,
    prep_iso_table(
      data = d,
      isotope_names = c("d13C", "d15N"),
      source_names = "Beaver",
      means_sds = FALSE
    )
  )
})

test_that("prep_iso_table (raw) errors on missing values", {
  d <- data.frame(Source = c("Beaver", "Beaver"), d13C = c(1, NA))
  expect_snapshot(
    error = TRUE,
    prep_iso_table(
      data = d,
      isotope_names = "d13C",
      source_names = "Beaver",
      means_sds = FALSE
    )
  )
})

# prep_iso_table(): summary mode ---------------------------------------------

test_that("prep_iso_table (summary) reorders rows to match source_names", {
  d <- data.frame(
    Source = c("Deer", "Beaver"),
    d13C_mean = c(-28, -26),
    d13C_sd = c(0.8, 1.0),
    d15N_mean = c(2, 3),
    d15N_sd = c(0.5, 0.6)
  )
  out <- prep_iso_table(
    data = d,
    isotope_names = c("d13C", "d15N"),
    source_names = c("Beaver", "Deer"),
    means_sds = TRUE
  )
  expect_equal(out$mode, "summary")
  expect_equal(unname(out$mean), matrix(c(-26, -28, 3, 2), nrow = 2))
  expect_equal(unname(out$sd), matrix(c(1.0, 0.8, 0.6, 0.5), nrow = 2))
})

test_that("prep_iso_table (summary) errors on missing mean/sd columns", {
  d <- data.frame(Source = "Beaver", d13C_mean = -26)
  expect_snapshot(
    error = TRUE,
    prep_iso_table(
      data = d,
      isotope_names = "d13C",
      source_names = "Beaver",
      means_sds = TRUE
    )
  )
})

test_that("prep_iso_table (summary) errors on duplicated source rows", {
  d <- data.frame(
    Source = c("Beaver", "Beaver"),
    d13C_mean = c(-26, -27),
    d13C_sd = c(1, 1)
  )
  expect_snapshot(
    error = TRUE,
    prep_iso_table(
      data = d,
      isotope_names = "d13C",
      source_names = "Beaver",
      means_sds = TRUE
    )
  )
})

test_that("prep_iso_table (summary) errors when source set doesn't match exactly", {
  d <- data.frame(Source = "Beaver", d13C_mean = -26, d13C_sd = 1)
  expect_snapshot(
    error = TRUE,
    prep_iso_table(
      data = d,
      isotope_names = "d13C",
      source_names = c("Beaver", "Deer"),
      means_sds = TRUE
    )
  )
})

test_that("prep_iso_table (summary) errors on missing values", {
  d <- data.frame(Source = "Beaver", d13C_mean = NA_real_, d13C_sd = 1)
  expect_snapshot(
    error = TRUE,
    prep_iso_table(
      data = d,
      isotope_names = "d13C",
      source_names = "Beaver",
      means_sds = TRUE
    )
  )
})

test_that("prep_iso_table (summary) errors on negative sd", {
  d <- data.frame(Source = "Beaver", d13C_mean = -26, d13C_sd = -1)
  expect_snapshot(
    error = TRUE,
    prep_iso_table(
      data = d,
      isotope_names = "d13C",
      source_names = "Beaver",
      means_sds = TRUE
    )
  )
})

# prep_conc_dep() -------------------------------------------------------------

test_that("prep_conc_dep returns NULL when conc_dep is FALSE", {
  expect_null(prep_conc_dep(
    source_data = data.frame(Source = "Beaver"),
    isotope_names = "d13C",
    source_names = "Beaver",
    source_means_sds = FALSE,
    conc_dep = FALSE
  ))
})

test_that("prep_conc_dep errors when conc_dep is not TRUE/FALSE", {
  expect_snapshot(
    error = TRUE,
    prep_conc_dep(
      source_data = data.frame(Source = "Beaver"),
      isotope_names = "d13C",
      source_names = "Beaver",
      source_means_sds = FALSE,
      conc_dep = "yes"
    )
  )
})

test_that("prep_conc_dep (raw) averages per-sample concentrations by source", {
  d <- data.frame(
    Source = c("Beaver", "Beaver", "Deer"),
    d13C_conc = c(0.3, 0.5, 0.4)
  )
  m <- prep_conc_dep(
    source_data = d,
    isotope_names = "d13C",
    source_names = c("Beaver", "Deer"),
    source_means_sds = FALSE,
    conc_dep = TRUE
  )
  expect_equal(dim(m), c(2L, 1L))
  expect_equal(unname(m[, 1]), c(0.4, 0.4))
})

test_that("prep_conc_dep (raw) handles multiple isotopes with the right shape", {
  d <- data.frame(
    Source = c("Beaver", "Beaver", "Deer"),
    d13C_conc = c(0.3, 0.5, 0.4),
    d15N_conc = c(0.05, 0.05, 0.04)
  )
  m <- prep_conc_dep(
    source_data = d,
    isotope_names = c("d13C", "d15N"),
    source_names = c("Beaver", "Deer"),
    source_means_sds = FALSE,
    conc_dep = TRUE
  )
  expect_equal(dim(m), c(2L, 2L))
  expect_equal(unname(m), matrix(c(0.4, 0.4, 0.05, 0.04), nrow = 2))
})

test_that("prep_conc_dep (summary) reads one concentration value per source", {
  d <- data.frame(Source = c("Deer", "Beaver"), d13C_conc = c(0.4, 0.35))
  m <- prep_conc_dep(
    source_data = d,
    isotope_names = "d13C",
    source_names = c("Beaver", "Deer"),
    source_means_sds = TRUE,
    conc_dep = TRUE
  )
  expect_equal(dim(m), c(2L, 1L))
  expect_equal(unname(m[, 1]), c(0.35, 0.4))
})

test_that("prep_conc_dep errors on missing concentration columns", {
  d <- data.frame(Source = "Beaver")
  expect_snapshot(
    error = TRUE,
    prep_conc_dep(
      source_data = d,
      isotope_names = "d13C",
      source_names = "Beaver",
      source_means_sds = FALSE,
      conc_dep = TRUE
    )
  )
})

test_that("prep_conc_dep (summary) errors on duplicated source rows", {
  d <- data.frame(Source = c("Beaver", "Beaver"), d13C_conc = c(0.3, 0.4))
  expect_snapshot(
    error = TRUE,
    prep_conc_dep(
      source_data = d,
      isotope_names = "d13C",
      source_names = "Beaver",
      source_means_sds = TRUE,
      conc_dep = TRUE
    )
  )
})

test_that("prep_conc_dep errors on missing concentration values", {
  d <- data.frame(Source = c("Beaver", "Beaver"), d13C_conc = c(0.3, NA))
  expect_snapshot(
    error = TRUE,
    prep_conc_dep(
      source_data = d,
      isotope_names = "d13C",
      source_names = "Beaver",
      source_means_sds = FALSE,
      conc_dep = TRUE
    )
  )
})

test_that("prep_conc_dep errors on values outside (0, 1]", {
  d <- data.frame(Source = "Beaver", d13C_conc = 1.5)
  expect_snapshot(
    error = TRUE,
    prep_conc_dep(
      source_data = d,
      isotope_names = "d13C",
      source_names = "Beaver",
      source_means_sds = TRUE,
      conc_dep = TRUE
    )
  )
})

test_that("prep_conc_dep errors when a source's concentrations sum to more than 1", {
  d <- data.frame(Source = "Beaver", d13C_conc = 0.7, d15N_conc = 0.6)
  expect_snapshot(
    error = TRUE,
    prep_conc_dep(
      source_data = d,
      isotope_names = c("d13C", "d15N"),
      source_names = "Beaver",
      source_means_sds = TRUE,
      conc_dep = TRUE
    )
  )
})
