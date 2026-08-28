mixture_data <- data.frame(
  d13C = c(-20, -21, -19, -22, -20.5, -21.5, -20.2, -21.2),
  d15N = c(10, 11, 9, 12, 10.5, 11.5, 10.3, 11.3),
  elevation = c(100, 120, 110, 300, 320, 310, 150, 280),
  rainfall = c(50, 55, 52, 80, 85, 82, 60, 78),
  Region = factor(c("A", "A", "A", "B", "B", "B", "A", "B")),
  Sex = factor(c("F", "M", "F", "M", "F", "M", "F", "M")),
  Pack = factor(c("P1", "P1", "P2", "P2", "P1", "P2", "P1", "P2"))
)
source_data <- data.frame(
  Source = c("Beaver", "Deer"),
  d13C_mean = c(-25, -18), d13C_sd = c(1, 1),
  d15N_mean = c(5, 8), d15N_sd = c(1, 1)
)
tdf_data <- data.frame(
  Source = c("Beaver", "Deer"),
  d13C_mean = c(1, 1.2), d13C_sd = c(0.2, 0.3),
  d15N_mean = c(3, 3.1), d15N_sd = c(0.4, 0.5)
)

skip_if_not_installed("cmdstanr")
skip_if_not_installed("ggplot2")

grDevices::pdf(nullfile())

fit <- bsimm(
  formula = ~ elevation + rainfall + Region + Sex + (1 | Pack), mixture_data = mixture_data,
  source_data = source_data, tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), source_means_sds = TRUE,
  chains = 1, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0,
  show_messages = FALSE, show_exceptions = FALSE
)
fit_interaction <- bsimm(
  formula = ~ elevation * Region + (1 | Pack), mixture_data = mixture_data,
  source_data = source_data, tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), source_means_sds = TRUE,
  chains = 1, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0,
  show_messages = FALSE, show_exceptions = FALSE
)
fit_no_fixed <- bsimm(
  formula = ~ 1 + (1 | Pack), mixture_data = mixture_data, source_data = source_data,
  tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), source_means_sds = TRUE,
  chains = 1, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0,
  show_messages = FALSE, show_exceptions = FALSE
)

test_that("conditional_effects errors on a non-bsimms_fit object", {
  expect_snapshot(error = TRUE, conditional_effects(list()))
})

test_that("conditional_effects errors when the model has no fixed-effect covariates", {
  expect_snapshot(error = TRUE, conditional_effects(fit_no_fixed))
})

test_that("conditional_effects errors when effects names an unknown covariate", {
  expect_snapshot(error = TRUE, conditional_effects(fit, effects = "banana"))
})

test_that("conditional_effects computes one element per fixed-effect covariate by default", {
  ce <- conditional_effects(fit)
  expect_s3_class(ce, "bsimms_conditional_effects")
  expect_equal(names(ce), c("elevation", "rainfall", "Region", "Sex"))
})

test_that("conditional_effects's default includes a two-way interaction present in the model", {
  ce <- conditional_effects(fit_interaction)
  expect_equal(names(ce), c("elevation", "Region", "elevation:Region"))
  df <- ce[["elevation:Region"]]
  expect_equal(attr(df, "moderator"), "Region")
})

test_that("conditional_effects's default drops (rather than errors on) a three-way+ interaction", {
  fit_3way <- bsimm(
    formula = ~ elevation * Region * Sex + (1 | Pack), mixture_data = mixture_data,
    source_data = source_data, tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), source_means_sds = TRUE,
    chains = 1, iter_warmup = 200, iter_sampling = 100, seed = 1, refresh = 0,
    show_messages = FALSE, show_exceptions = FALSE
  )
  ce <- conditional_effects(fit_3way)
  expect_setequal(
    names(ce),
    c("elevation", "Region", "Sex", "elevation:Region", "elevation:Sex", "Region:Sex")
  )
})

test_that("conditional_effects respects an explicit effects subset", {
  ce <- conditional_effects(fit, effects = "elevation")
  expect_equal(names(ce), "elevation")
})

test_that("conditional_effects builds a resolution-point grid for a numeric covariate", {
  ce <- conditional_effects(fit, effects = "elevation", resolution = 10)
  expect_equal(length(unique(ce$elevation$row)), 10)
  expect_true(attr(ce$elevation, "is_numeric"))
  expect_equal(names(ce$elevation), c("row", "source", "estimate", "lower", "upper", "width", "elevation"))
})

test_that("conditional_effects uses observed levels for a factor covariate", {
  ce <- conditional_effects(fit, effects = "Region")
  expect_equal(length(unique(ce$Region$row)), 2)
  expect_false(attr(ce$Region, "is_numeric"))
  expect_setequal(unique(ce$Region$Region), c("A", "B"))
})

test_that("build_conditional_grid holds a factor's reference at its first level, not the mode", {
  fixed_frame <- data.frame(
    elevation = c(100, 200, 300),
    Region = factor(c("A", "A", "B"), levels = c("B", "A"))
  )
  cond <- build_conditional_grid(fixed_frame, "elevation", resolution = 5, int_conditions = NULL, ref_conditions = NULL)
  expect_equal(as.character(unique(cond$grid$Region)), "B")
})

test_that("ref_conditions overrides the default reference value for a numeric covariate", {
  ce <- conditional_effects(fit, effects = "Region", ref_conditions = list(elevation = 500))
  # can't read elevation off ce$Region directly (it's held constant, not attached to the
  # output), so check indirectly via build_conditional_grid instead
  cond <- build_conditional_grid(
    fit$spec$fixed_frame, "Region", resolution = 5,
    int_conditions = NULL, ref_conditions = list(elevation = 500)
  )
  expect_equal(unique(cond$grid$elevation), 500)
})

test_that("ref_conditions overrides the default reference level for a factor covariate", {
  cond <- build_conditional_grid(
    fit$spec$fixed_frame, "elevation", resolution = 5,
    int_conditions = NULL, ref_conditions = list(Region = "B")
  )
  expect_equal(as.character(unique(cond$grid$Region)), "B")
})

test_that("ref_conditions errors on non-numeric values for a numeric covariate", {
  expect_snapshot(
    error = TRUE,
    conditional_effects(fit, effects = "Region", ref_conditions = list(elevation = "high"))
  )
})

test_that("ref_conditions errors on an unseen level for a factor covariate", {
  expect_snapshot(
    error = TRUE,
    conditional_effects(fit, effects = "elevation", ref_conditions = list(Region = "X"))
  )
})

test_that("ref_conditions errors on more than one value for a covariate", {
  expect_snapshot(
    error = TRUE,
    conditional_effects(fit, effects = "Region", ref_conditions = list(elevation = c(100, 200)))
  )
})

test_that("ref_conditions errors on a misspelled covariate name instead of falling back to the default", {
  expect_snapshot(
    error = TRUE,
    conditional_effects(fit, effects = "Region", ref_conditions = list(elevaton = 500))
  )
})

test_that("int_conditions errors on a misspelled covariate name instead of falling back to the default", {
  expect_snapshot(
    error = TRUE,
    conditional_effects(fit, effects = "elevation:Region", int_conditions = list(Regoin = "A"))
  )
})

test_that("conditional_effects's ... is forwarded to posterior_proportions", {
  expect_snapshot(error = TRUE, conditional_effects(fit, effects = "elevation", ndraws = 1000))
})

test_that("conditional_effects errors when re_formula names a term absent from the grid", {
  expect_snapshot(error = TRUE, conditional_effects(fit, effects = "elevation", re_formula = ~ (1 | Pack)))
})

test_that("plot.bsimms_conditional_effects returns one ggplot per covariate", {
  ce <- conditional_effects(fit, effects = c("elevation", "Region"))
  plots <- plot(ce)
  expect_type(plots, "list")
  expect_equal(names(plots), c("elevation", "Region"))
  expect_s3_class(plots$elevation, "ggplot")
  expect_s3_class(plots$Region, "ggplot")
})

test_that("printing a bsimms_conditional_effects object draws its plots", {
  ce <- conditional_effects(fit, effects = "elevation")
  expect_invisible(print(ce))
})

test_that("plot's ... is forwarded to geom_ribbon for a numeric covariate", {
  ce <- conditional_effects(fit, effects = "elevation")
  plots <- plot(ce, linetype = "dashed")
  expect_equal(plots$elevation$layers[[1]]$aes_params$linetype, "dashed")
})

test_that("plot's ... is forwarded to geom_linerange for a factor covariate", {
  ce <- conditional_effects(fit, effects = "Region")
  plots <- plot(ce, linetype = "dashed")
  expect_equal(plots$Region$layers[[1]]$aes_params$linetype, "dashed")
})

test_that("plot's plot = FALSE builds the plots without displaying them", {
  ce <- conditional_effects(fit, effects = "elevation")
  plots <- plot(ce, plot = FALSE)
  expect_s3_class(plots$elevation, "ggplot")
})

test_that("a plot built with plot = FALSE can be further customised before printing", {
  ce <- conditional_effects(fit, effects = "elevation")
  plots <- plot(ce, plot = FALSE)
  customised <- plots$elevation + ggplot2::labs(title = "Custom title")
  expect_s3_class(customised, "ggplot")
  expect_equal(customised$labels$title, "Custom title")
})

# two-way interactions --------------------------------------------------------

test_that("conditional_effects errors when an interaction names more than two covariates", {
  expect_snapshot(error = TRUE, conditional_effects(fit, effects = "elevation:Region:Sex"))
})

test_that("conditional_effects errors when an interaction names the same covariate twice", {
  expect_snapshot(error = TRUE, conditional_effects(fit, effects = "elevation:elevation"))
})

test_that("conditional_effects errors when an interaction names an unknown covariate", {
  expect_snapshot(error = TRUE, conditional_effects(fit, effects = "elevation:banana"))
})

test_that("conditional_effects builds a continuous:factor interaction grid, faceted by the factor moderator", {
  ce <- conditional_effects(fit, effects = "elevation:Region", resolution = 5)
  df <- ce[["elevation:Region"]]
  expect_true(attr(df, "is_numeric"))
  expect_equal(attr(df, "moderator"), "Region")
  expect_false(attr(df, "moderator_is_numeric"))
  expect_setequal(unique(df$Region), c("A", "B"))
  expect_equal(length(unique(df$elevation)), 5)
})

test_that("conditional_effects builds a factor:factor interaction grid", {
  ce <- conditional_effects(fit, effects = "Region:Sex")
  df <- ce[["Region:Sex"]]
  expect_false(attr(df, "is_numeric"))
  expect_equal(attr(df, "moderator"), "Sex")
  expect_setequal(unique(df$Region), c("A", "B"))
  expect_setequal(unique(df$Sex), c("F", "M"))
})

test_that("conditional_effects builds a continuous:continuous interaction with mean +/- 1 SD by default", {
  ce <- conditional_effects(fit, effects = "elevation:rainfall", resolution = 5)
  df <- ce[["elevation:rainfall"]]
  expect_true(attr(df, "moderator_is_numeric"))
  expect_equal(length(unique(df$rainfall)), 3)
  expect_equal(
    sort(unique(df$rainfall)),
    sort(mean(mixture_data$rainfall) + c(-1, 0, 1) * sd(mixture_data$rainfall))
  )
})

test_that("int_conditions overrides the default representative values for a numeric moderator", {
  ce <- conditional_effects(
    fit,
    effects = "elevation:rainfall", resolution = 5, int_conditions = list(rainfall = c(50, 80))
  )
  df <- ce[["elevation:rainfall"]]
  expect_setequal(unique(df$rainfall), c(50, 80))
})

test_that("int_conditions overrides the default levels for a factor moderator", {
  ce <- conditional_effects(fit, effects = "Region:Sex", int_conditions = list(Sex = "F"))
  df <- ce[["Region:Sex"]]
  expect_setequal(unique(df$Sex), "F")
})

test_that("int_conditions errors on an unseen level for a factor moderator", {
  expect_snapshot(error = TRUE, conditional_effects(fit, effects = "Region:Sex", int_conditions = list(Sex = "X")))
})

test_that("int_conditions errors on non-numeric values for a numeric moderator", {
  expect_snapshot(
    error = TRUE,
    conditional_effects(fit, effects = "elevation:rainfall", int_conditions = list(rainfall = "high"))
  )
})

test_that("plot facets by the moderator for a two-way interaction", {
  ce <- conditional_effects(fit, effects = "elevation:Region", resolution = 5)
  plots <- plot(ce)
  built <- ggplot2::ggplot_build(plots[["elevation:Region"]])
  expect_equal(length(unique(built$data[[1]]$PANEL)), 2)
})

# method = posterior_epred / posterior_predict -------------------------------

test_that("conditional_effects errors on an invalid method", {
  expect_snapshot(error = TRUE, conditional_effects(fit, effects = "elevation", method = "banana"))
})

test_that("conditional_effects requires resp for a multi-isotope model with method = posterior_epred", {
  skip_if_not_installed("rstantools")
  expect_snapshot(error = TRUE, conditional_effects(fit, effects = "elevation", method = "posterior_epred"))
})

test_that("conditional_effects errors on an unknown resp", {
  skip_if_not_installed("rstantools")
  expect_snapshot(
    error = TRUE,
    conditional_effects(fit, effects = "elevation", method = "posterior_epred", resp = "banana")
  )
})

test_that("method = posterior_epred computes an isotope-labelled conditional effect", {
  skip_if_not_installed("rstantools")
  ce <- conditional_effects(fit, effects = "elevation", method = "posterior_epred", resp = "d13C")
  df <- ce$elevation
  expect_equal(attr(df, "cat_col"), "isotope")
  expect_equal(attr(df, "y_label"), "d13C")
  expect_equal(unique(df$isotope), "d13C")
})

test_that("method = posterior_predict has wider intervals than posterior_epred", {
  skip_if_not_installed("rstantools")
  ce_epred <- conditional_effects(
    fit,
    effects = "elevation", method = "posterior_epred", resp = "d13C", probs = 0.95
  )
  ce_predict <- conditional_effects(
    fit,
    effects = "elevation", method = "posterior_predict", resp = "d13C", probs = 0.95
  )
  width_epred <- mean(ce_epred$elevation$upper - ce_epred$elevation$lower)
  width_predict <- mean(ce_predict$elevation$upper - ce_predict$elevation$lower)
  expect_gt(width_predict, width_epred)
})

test_that("plot has no colour legend for method = posterior_epred/posterior_predict (single isotope)", {
  skip_if_not_installed("rstantools")
  ce <- conditional_effects(fit, effects = "elevation", method = "posterior_epred", resp = "d13C")
  g <- plot(ce, plot = FALSE)$elevation
  expect_null(g$labels$colour)
})

test_that("plot still has a colour legend for the default method = posterior_proportions", {
  ce <- conditional_effects(fit, effects = "elevation")
  g <- plot(ce, plot = FALSE)$elevation
  expect_equal(g$labels$colour, "Source")
})

test_that("interval dodging doesn't split nested widths apart when there's no legend", {
  skip_if_not_installed("rstantools")
  ce <- conditional_effects(fit, effects = "Region", method = "posterior_epred", resp = "d13C")
  g <- plot(ce, plot = FALSE)$Region
  built <- ggplot2::ggplot_build(g)
  expect_equal(length(unique(built$data[[1]]$x)), 2)
})
