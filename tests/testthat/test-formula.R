test_that("parse_bsimms_formula extracts fixed and random terms", {
  d <- data.frame(
    Sex = factor(c("M", "F", "F", "M", "M", "F")),
    Region = factor(c("A", "A", "B", "B", "C", "C")),
    x = c(1, 2, 3, 4, 5, 6)
  )
  pf <- parse_bsimms_formula(formula = ~ Sex + x + (1 | Region), data = d)
  expect_true("(Intercept)" %in% pf$fixed_names)
  expect_true(any(grepl("^Sex", pf$fixed_names)))
  expect_equal(nrow(pf$X), 6)
  expect_length(pf$re_terms, 1)
  expect_equal(pf$re_terms[[1]]$group, "Region")
  expect_equal(length(pf$re_terms[[1]]$group_levels), 3)
})

test_that("a character fixed-effect column is coerced to factor", {
  d <- data.frame(Sex = c("M", "F", "F", "M", "M", "F"))
  pf <- parse_bsimms_formula(formula = ~Sex, data = d)
  expect_true(is.factor(pf$fixed_frame$Sex))
})

test_that("double-bar terms expand into uncorrelated single-bar terms", {
  d <- data.frame(Region = factor(rep(c("A", "B"), each = 5)), x = rnorm(10))
  pf <- parse_bsimms_formula(formula = ~ x + (x || Region), data = d)
  expect_length(pf$re_terms, 2)
})

test_that("a two-sided formula is rejected", {
  d <- data.frame(y = 1:3, x = 1:3)
  expect_snapshot(error = TRUE, parse_bsimms_formula(formula = y ~ x, data = d))
})

test_that("a non-formula is rejected", {
  d <- data.frame(y = 1:3, x = 1:3)
  expect_snapshot(error = TRUE, parse_bsimms_formula(formula = "~x", data = d))
})

test_that("non-data-frame data is rejected", {
  expect_snapshot(
    error = TRUE,
    parse_bsimms_formula(formula = ~x, data = list(x = 1:3))
  )
})

test_that("a chained bar term is rejected", {
  d <- data.frame(
    a = factor(rep(c("x", "y"), 5)),
    b = factor(rep(c("p", "q"), 5))
  )
  expect_snapshot(
    error = TRUE,
    parse_bsimms_formula(formula = ~ 1 + (1 | a | b), data = d)
  )
})

test_that("a grouping factor with missing values is rejected", {
  d <- data.frame(Region = factor(c("A", "A", "B", NA)), x = c(1, 2, 3, 4))
  expect_snapshot(
    error = TRUE,
    parse_bsimms_formula(formula = ~ 1 + (1 | Region), data = d)
  )
})

test_that("fixed-effect interaction terms are captured correctly", {
  d <- data.frame(
    Sex = factor(c("M", "F", "F", "M", "M", "F")),
    Region = factor(c("A", "A", "B", "B", "C", "C"))
  )
  pf <- parse_bsimms_formula(formula = ~ Sex * Region, data = d)
  expect_true(any(grepl("^SexM:Region", pf$fixed_names)))
  expect_equal(ncol(pf$X), ncol(stats::model.matrix(~ Sex * Region, data = d)))
})

test_that("crossed group-level terms (a:b) are captured correctly", {
  # a, b left as character (not factor) to also exercise the fact that `:`
  # only computes the factor interaction when both operands are factors.
  d <- data.frame(a = rep(c("x", "y"), each = 4), b = rep(c("p", "q"), 4))
  pf <- parse_bsimms_formula(formula = ~ 1 + (1 | a:b), data = d)
  expect_equal(pf$re_terms[[1]]$group, "a:b")
  expect_setequal(pf$re_terms[[1]]$group_levels, c("x:p", "x:q", "y:p", "y:q"))
})

test_that("nested group-level terms (site + site:individual) are captured correctly", {
  # individual labels are reused across sites, so site:individual (not
  # individual alone) is needed to identify each nested unit uniquely.
  d <- data.frame(
    site = rep(c("A", "B"), each = 4),
    individual = rep(rep(c("1", "2"), each = 2), 2)
  )
  pf <- parse_bsimms_formula(
    formula = ~ 1 + (1 | site) + (1 | site:individual),
    data = d
  )
  expect_length(pf$re_terms, 2)
  expect_equal(pf$re_terms[[1]]$group, "site")
  expect_length(pf$re_terms[[1]]$group_levels, 2)
  expect_equal(pf$re_terms[[2]]$group, "site:individual")
  expect_setequal(pf$re_terms[[2]]$group_levels, c("A:1", "A:2", "B:1", "B:2"))
})

test_that("the site/individual shorthand is expanded the same way as lme4", {
  d <- data.frame(
    site = rep(c("A", "B"), each = 4),
    individual = rep(rep(c("1", "2"), each = 2), 2)
  )
  pf <- parse_bsimms_formula(formula = ~ 1 + (1 | site / individual), data = d)
  expect_length(pf$re_terms, 2)
  groups <- vapply(pf$re_terms, function(x) x$group, character(1))
  expect_setequal(groups, c("individual:site", "site"))
})

test_that("continuous columns are classified correctly for standardisation", {
  d <- data.frame(
    Age = c(1, 2, 3, 4, 5, 6),
    Weight = c(10, 20, 30, 40, 50, 60),
    Sex = factor(c("M", "F", "M", "F", "M", "F")),
    Flag = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
  )
  classify <- function(f) {
    pf <- parse_bsimms_formula(formula = f, data = d)
    stats::setNames(pf$x_continuous, pf$fixed_names)
  }

  expect_equal(classify(~1), c("(Intercept)" = FALSE))
  expect_equal(classify(~Age)[["Age"]], TRUE)
  expect_false(any(classify(~Sex)[-1]))
  expect_true(classify(~ Age:Weight)[["Age:Weight"]])
  expect_false(any(classify(~ Age:Sex)[-1]))
  expect_false(any(classify(~Flag)[-1]))

  mixed <- classify(~ Age + Sex + Age:Sex)
  expect_true(mixed[["Age"]])
  expect_false(mixed[[grep("^Sex", names(mixed), value = TRUE)]])
  expect_false(mixed[[grep("^Age:Sex", names(mixed), value = TRUE)]])
})
