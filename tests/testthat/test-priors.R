# bsimms_prior() --------------------------------------------------------

test_that("bsimms_prior builds a one-row bsimms_prior data frame", {
  p <- bsimms_prior(prior = "normal(0, 2)", class = "b")
  expect_s3_class(p, "bsimms_prior")
  expect_equal(p$prior, "normal(0, 2)")
  expect_equal(p$class, "b")
  expect_equal(p$coef, "")
  expect_equal(p$resp, "")
  expect_equal(p$group, "")
})

test_that("bsimms_prior errors on an invalid class", {
  expect_snapshot(error = TRUE, bsimms_prior(prior = "normal(0, 1)", class = "banana"))
})

test_that("bsimms_prior errors on a non-string prior", {
  expect_snapshot(error = TRUE, bsimms_prior(prior = 1, class = "b"))
  expect_snapshot(error = TRUE, bsimms_prior(prior = c("a", "b"), class = "b"))
})

# c.bsimms_prior() --------------------------------------------------------

test_that("c.bsimms_prior row-binds multiple priors", {
  p <- c(
    bsimms_prior(prior = "normal(0, 2)", class = "b"),
    bsimms_prior(prior = "normal(1, 1)", class = "b", coef = "SeasonWinter")
  )
  expect_s3_class(p, "bsimms_prior")
  expect_equal(nrow(p), 2)
  expect_equal(p$coef, c("", "SeasonWinter"))
})

# print.bsimms_prior() --------------------------------------------------------

test_that("print.bsimms_prior prints a table with headers", {
  p <- c(
    bsimms_prior(prior = "normal(0, 2)", class = "b"),
    bsimms_prior(prior = "student_t(3, 0, 1)", class = "sd", group = "Region")
  )
  expect_snapshot(print(p))
})

# merge_bsimms_prior() --------------------------------------------------------

test_that("merge_bsimms_prior overrides matching rows and appends new ones", {
  default <- c(
    bsimms_prior(prior = "normal(0, 1)", class = "b"),
    bsimms_prior(prior = "student_t(3, 0, 2)", class = "sd", group = "Region")
  )
  user <- c(
    bsimms_prior(prior = "normal(0, 5)", class = "b"),
    bsimms_prior(prior = "normal(3.4, 0.3)", class = "tdf_mean", resp = "d15N", group = "Beaver")
  )
  merged <- merge_bsimms_prior(default, user)
  expect_equal(nrow(merged), 3)
  expect_equal(merged$prior[merged$class == "b"], "normal(0, 5)")
  expect_equal(merged$prior[merged$class == "tdf_mean"], "normal(3.4, 0.3)")
})

test_that("merge_bsimms_prior returns default unchanged when user is NULL", {
  default <- bsimms_prior(prior = "normal(0, 1)", class = "b")
  expect_identical(merge_bsimms_prior(default, NULL), default)
})

test_that("merge_bsimms_prior errors on a non-bsimms_prior user argument", {
  default <- bsimms_prior(prior = "normal(0, 1)", class = "b")
  expect_snapshot(error = TRUE, merge_bsimms_prior(default, data.frame(x = 1)))
})

# select_prior() --------------------------------------------------------

test_that("select_prior falls back from most to least specific", {
  df <- c(
    bsimms_prior(prior = "student_t(3, 0, 1)", class = "sd"),
    bsimms_prior(prior = "student_t(3, 0, 2)", class = "sd", group = "Region"),
    bsimms_prior(prior = "student_t(3, 0, 3)", class = "sd", coef = "x", group = "Region")
  )
  expect_equal(select_prior(df, "sd", coef = "x", group = "Region"), "student_t(3, 0, 3)")
  expect_equal(select_prior(df, "sd", coef = "y", group = "Region"), "student_t(3, 0, 2)")
  expect_equal(select_prior(df, "sd", coef = "y", group = "Other"), "student_t(3, 0, 1)")
})

test_that("select_prior falls back across coef, resp and group jointly", {
  df <- c(
    bsimms_prior(prior = "student_t(3, 0, 1)", class = "source_sd"),
    bsimms_prior(prior = "student_t(3, 0, 2)", class = "source_sd", group = "Beaver"),
    bsimms_prior(prior = "student_t(3, 0, 3)", class = "source_sd", resp = "d13C"),
    bsimms_prior(prior = "student_t(3, 0, 4)", class = "source_sd", resp = "d13C", group = "Beaver")
  )
  expect_equal(select_prior(df, "source_sd", resp = "d13C", group = "Beaver"), "student_t(3, 0, 4)")
  expect_equal(select_prior(df, "source_sd", resp = "d13C", group = "Deer"), "student_t(3, 0, 3)")
  expect_equal(select_prior(df, "source_sd", resp = "d15N", group = "Beaver"), "student_t(3, 0, 2)")
  expect_equal(select_prior(df, "source_sd", resp = "d15N", group = "Deer"), "student_t(3, 0, 1)")
})

test_that("select_prior errors when the class is entirely absent", {
  df <- bsimms_prior(prior = "normal(0, 1)", class = "b")
  expect_snapshot(error = TRUE, select_prior(df, "sd"))
})

test_that("select_prior errors when no row matches any specificity level", {
  df <- bsimms_prior(prior = "normal(0, 1)", class = "sd", coef = "x", group = "A")
  expect_snapshot(error = TRUE, select_prior(df, "sd", coef = "y", group = "B"))
})

# default_bsimms_prior() / bsimms_get_prior() ----------------------------

mixture_data <- data.frame(
  d13C = c(-20, -21, -19, -22, -20.5, -21.5),
  d15N = c(10, 11, 9, 12, 10.5, 11.5),
  Region = factor(rep(c("A", "B"), each = 3))
)
source_data <- data.frame(
  Source = c("Beaver", "Beaver", "Beaver", "Beaver", "Beaver", "Deer", "Deer", "Deer"),
  d13C = c(-26, -25, -24, -25, -10, -18, -17, -19),
  d15N = c(5, 6, 4, 5.5, 5, 8, 9, 7)
)
tdf_data <- data.frame(
  Source = c("Beaver", "Deer"),
  d13C_mean = c(1, 1.2), d13C_sd = c(0.2, 0.3),
  d15N_mean = c(3, 3.1), d15N_sd = c(0.4, 0.5)
)

test_that("bsimms_get_prior returns default priors for a simple model", {
  p <- bsimms_get_prior(
    formula = ~1, mixture_data = mixture_data, source_data = source_data,
    tdf_data = tdf_data, isotope_names = c("d13C", "d15N")
  )
  expect_s3_class(p, "bsimms_prior")
  expect_true(all(c("b", "p_global", "resid_prop", "source_mean", "source_sd") %in% p$class))
  expect_equal(sort(p$group[p$class == "p_global"]), c("Beaver", "Deer"))
})

test_that("default source_mean/source_sd priors are scaled by median/MAD, not mean/SD", {
  p <- bsimms_get_prior(
    formula = ~1, mixture_data = mixture_data, source_data = source_data,
    tdf_data = tdf_data, isotope_names = c("d13C", "d15N")
  )
  med <- stats::median(source_data$d13C)
  mad <- stats::mad(source_data$d13C)
  expect_false(isTRUE(all.equal(med, mean(source_data$d13C))))
  expect_equal(
    p$prior[p$class == "source_mean" & p$resp == "d13C"],
    sprintf("normal(%.6g, %.6g)", med, 10 * mad)
  )
  expect_equal(
    p$prior[p$class == "source_sd" & p$resp == "d13C"],
    sprintf("student_t(3, 0, %.6g)", mad)
  )
})

test_that("default sigma prior (residual_only) is scaled by median/MAD", {
  p <- bsimms_get_prior(
    formula = ~1, mixture_data = mixture_data, source_data = source_data,
    tdf_data = tdf_data, isotope_names = c("d13C", "d15N"),
    error_structure = "residual_only"
  )
  d13c_mad <- stats::mad(mixture_data$d13C)
  expect_equal(
    p$prior[p$class == "sigma" & p$resp == "d13C"],
    sprintf("student_t(3, 0, %.6g)", d13c_mad)
  )
})
