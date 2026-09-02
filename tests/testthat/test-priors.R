mixture_data <- data.frame(
  d13C = c(-20, -21, -19, -22, -20.5, -21.5),
  d15N = c(10, 11, 9, 12, 10.5, 11.5),
  Region = factor(rep(c("A", "B"), each = 3)),
  Sex = factor(rep(c("M", "F"), 3))
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

spec <- build_bsimms_spec(
  formula = ~1, mixture_data = mixture_data, source_data = source_data,
  tdf_data = tdf_data, isotope_names = c("d13C", "d15N")
)
spec_re <- build_bsimms_spec(
  formula = ~ Sex + (1 | Region), mixture_data = mixture_data, source_data = source_data,
  tdf_data = tdf_data, isotope_names = c("d13C", "d15N")
)

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
  merged <- merge_bsimms_prior(default, user, spec_re)
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

test_that("a class-level override cascades onto untouched, more specific default rows", {
  default <- default_bsimms_prior(spec)
  user <- bsimms_prior("uniform(0, 5)", class = "resid_prop")
  merged <- merge_bsimms_prior(default, user, spec)
  rows <- merged[merged$class == "resid_prop", ]
  expect_true(all(rows$prior == "uniform(0, 5)"))
  expect_setequal(rows$resp, c("", "d13C", "d15N"))
})

test_that("a more specific override wins over a cascaded general one, regardless of order", {
  default <- default_bsimms_prior(spec)
  user_general_first <- c(
    bsimms_prior("uniform(0, 5)", class = "resid_prop"),
    bsimms_prior("uniform(0, 3)", class = "resid_prop", resp = "d13C")
  )
  user_specific_first <- c(
    bsimms_prior("uniform(0, 3)", class = "resid_prop", resp = "d13C"),
    bsimms_prior("uniform(0, 5)", class = "resid_prop")
  )
  merged1 <- merge_bsimms_prior(default, user_general_first, spec)
  merged2 <- merge_bsimms_prior(default, user_specific_first, spec)
  expect_identical(merged1$prior[merged1$class == "resid_prop"], merged2$prior[merged2$class == "resid_prop"])
  rows <- merged1[merged1$class == "resid_prop", ]
  expect_equal(rows$prior[rows$resp == "d13C"], "uniform(0, 3)")
  expect_equal(rows$prior[rows$resp == "d15N"], "uniform(0, 5)")
})

test_that("a group-only override cascades across every source for p_global", {
  default <- default_bsimms_prior(spec)
  user <- bsimms_prior("3", class = "p_global")
  merged <- merge_bsimms_prior(default, user, spec)
  rows <- merged[merged$class == "p_global" & merged$group != "", ]
  expect_true(all(rows$prior == "3"))
})

test_that("a resp-only override cascades onto every source, leaving other isotopes alone", {
  default <- default_bsimms_prior(spec)
  user <- bsimms_prior("normal(0, 50)", class = "source_mean", resp = "d13C")
  merged <- merge_bsimms_prior(default, user, spec)
  rows <- merged[merged$class == "source_mean" & merged$resp == "d13C" & merged$group != "", ]
  expect_true(all(rows$prior == "normal(0, 50)"))
  untouched <- default[default$class == "source_mean" & default$resp == "d15N", ]
  merged_d15n <- merged[merged$class == "source_mean" & merged$resp == "d15N" & merged$group != "", ]
  expect_equal(merged_d15n$prior, untouched$prior)
})

test_that("merge_bsimms_prior errors on an unrecognised coef/resp/group", {
  default <- default_bsimms_prior(spec_re)
  expect_snapshot(
    error = TRUE,
    merge_bsimms_prior(default, bsimms_prior("normal(0, 1)", class = "b", coef = "SexZZZ"), spec_re)
  )
  expect_snapshot(
    error = TRUE,
    merge_bsimms_prior(default, bsimms_prior("2", class = "p_global", group = "Beavr"), spec_re)
  )
  expect_snapshot(
    error = TRUE,
    merge_bsimms_prior(default, bsimms_prior("student_t(3, 0, 1)", class = "sd", group = "Regionn"), spec_re)
  )
})

test_that("merge_bsimms_prior errors clearly when a class has no valid terms in this model", {
  default <- default_bsimms_prior(spec)
  expect_snapshot(
    error = TRUE,
    merge_bsimms_prior(default, bsimms_prior("student_t(3, 0, 1)", class = "sd", group = "Region"), spec)
  )
})

test_that("merge_bsimms_prior accepts valid coef/resp/group without error", {
  default <- default_bsimms_prior(spec_re)
  merged <- merge_bsimms_prior(default, bsimms_prior("normal(0, 1)", class = "b", coef = "SexM"), spec_re)
  expect_equal(merged$prior[merged$class == "b" & merged$coef == "SexM"], "normal(0, 1)")
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

test_that("bsimms_get_prior returns default priors for a simple model", {
  p <- bsimms_get_prior(
    formula = ~1, mixture_data = mixture_data, source_data = source_data,
    tdf_data = tdf_data, isotope_names = c("d13C", "d15N")
  )
  expect_s3_class(p, "bsimms_prior")
  expect_true(all(c("b", "p_global", "resid_prop", "source_mean", "source_sd") %in% p$class))
  expect_equal(sort(p$group[p$class == "p_global"]), c("Beaver", "Deer"))
})

test_that("default source_mean/source_sd priors use each source's own median/MAD, not mean/SD", {
  p <- bsimms_get_prior(
    formula = ~1, mixture_data = mixture_data, source_data = source_data,
    tdf_data = tdf_data, isotope_names = c("d13C", "d15N")
  )
  beaver_d13c <- source_data$d13C[source_data$Source == "Beaver"]
  med <- stats::median(beaver_d13c)
  mad <- stats::mad(beaver_d13c)
  expect_false(isTRUE(all.equal(med, mean(beaver_d13c))))
  expect_equal(
    p$prior[p$class == "source_mean" & p$resp == "d13C" & p$group == "Beaver"],
    sprintf("normal(%.6g, %.6g)", med, 10 * mad)
  )
  expect_equal(
    p$prior[p$class == "source_sd" & p$resp == "d13C" & p$group == "Beaver"],
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
