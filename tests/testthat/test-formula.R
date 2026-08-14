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

test_that("double-bar terms expand into uncorrelated single-bar terms", {
  d <- data.frame(Region = factor(rep(c("A", "B"), each = 5)), x = rnorm(10))
  pf <- parse_bsimms_formula(formula = ~ x + (x || Region), data = d)
  expect_length(pf$re_terms, 2)
})

test_that("a two-sided formula is rejected", {
  d <- data.frame(y = 1:3, x = 1:3)
  expect_error(parse_bsimms_formula(formula = y ~ x, data = d), "left-hand side")
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
  pf <- parse_bsimms_formula(formula = ~ 1 + (1 | site) + (1 | site:individual), data = d)
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
