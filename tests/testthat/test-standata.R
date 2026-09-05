mixture_data <- data.frame(
  d13C = c(-20, -21, -19, -22, -20.5, -21.5),
  d15N = c(10, 11, 9, 12, 10.5, 11.5),
  Sex = factor(rep(c("M", "F"), 3)),
  Region = factor(rep(c("A", "B"), each = 3))
)
source_data_raw <- data.frame(
  Source = rep(c("Beaver", "Deer"), each = 3),
  d13C = c(-25, -24, -26, -18, -17, -19),
  d15N = c(5, 6, 4, 8, 9, 7)
)
source_data_summary <- data.frame(
  Source = c("Beaver", "Deer"),
  d13C_mean = c(-25, -18),
  d13C_sd = c(1, 1),
  d15N_mean = c(5, 8),
  d15N_sd = c(1, 1)
)
tdf_data_summary <- data.frame(
  Source = c("Beaver", "Deer"),
  d13C_mean = c(1, 1.2),
  d13C_sd = c(0.2, 0.3),
  d15N_mean = c(3, 3.1),
  d15N_sd = c(0.4, 0.5)
)

test_that("make_standata returns a bsimms_standata list", {
  sdata <- make_standata(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )
  expect_s3_class(sdata, "bsimms_standata")
  expect_type(sdata, "list")
  expect_equal(sdata$N, 6)
  expect_equal(sdata$J, 2)
  expect_equal(sdata$K, 2)
  expect_equal(sdata$D, 1)
  expect_equal(sdata$alpha_dirichlet, c(1, 1))
})

test_that("P/X are included only when the model has fixed effects", {
  sdata_no_fe <- make_standata(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )
  expect_null(sdata_no_fe$P)
  expect_null(sdata_no_fe$X)

  sdata_fe <- make_standata(
    formula = ~Sex,
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )
  expect_equal(sdata_fe$P, 1)
  expect_equal(dim(sdata_fe$X), c(6, 1))
})

test_that("group-level terms add N_re_*/M_re_*/Z_re_*/grp_re_* entries", {
  sdata <- make_standata(
    formula = ~ 1 + (1 | Region),
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )
  expect_equal(sdata$N_re_Region, 2)
  expect_equal(sdata$M_re_Region, 1)
  expect_equal(dim(sdata$Z_re_Region), c(6, 1))
  expect_equal(sdata$grp_re_Region, c(1, 1, 1, 2, 2, 2))
})

test_that("raw source data adds N_source_raw/source_idx/source_raw", {
  sdata <- make_standata(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )
  expect_equal(sdata$N_source_raw, 6)
  expect_equal(sdata$source_idx, c(1, 1, 1, 2, 2, 2))
  expect_equal(dim(sdata$source_raw), c(6, 2))
  expect_null(sdata$source_mean_data)
})

test_that("summarised source data adds source_mean_data/source_sd_data instead", {
  sdata <- make_standata(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_summary,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N"),
    source_means_sds = TRUE
  )
  expect_equal(dim(sdata$source_mean_data), c(2, 2))
  expect_equal(dim(sdata$source_sd_data), c(2, 2))
  expect_null(sdata$source_raw)
})

test_that("conc_dep = TRUE adds the conc matrix", {
  d <- source_data_raw
  d$d13C_conc <- c(0.3, 0.32, 0.28, 0.35, 0.33, 0.37)
  d$d15N_conc <- c(0.05, 0.05, 0.06, 0.04, 0.05, 0.04)
  sdata <- make_standata(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = d,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N"),
    conc_dep = TRUE
  )
  expect_equal(dim(sdata$conc), c(2, 2))
})

test_that("a user prior override changes alpha_dirichlet", {
  sdata_default <- make_standata(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )
  expect_equal(sdata_default$alpha_dirichlet, c(1, 1))

  sdata_override <- make_standata(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N"),
    prior = bsimms_prior(prior = "3", class = "p_global", group = "Beaver")
  )
  expect_equal(sdata_override$alpha_dirichlet, c(3, 1))
})

test_that("make_standata's data list is compatible with make_stancode's generated program", {
  skip_if_not_installed("cmdstanr")
  code <- make_stancode(
    formula = ~ Sex + (1 | Region),
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )
  sdata <- make_standata(
    formula = ~ Sex + (1 | Region),
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )
  f <- tempfile(fileext = ".stan")
  on.exit(unlink(f))
  writeLines(as.character(code), f)
  mod <- cmdstanr::cmdstan_model(f, compile = TRUE, quiet = TRUE)
  fit <- mod$sample(
    data = unclass(sdata),
    chains = 1,
    iter_warmup = 100,
    iter_sampling = 50,
    refresh = 0,
    seed = 1,
    show_messages = FALSE,
    show_exceptions = FALSE
  )
  draws <- fit$draws("p_global")
  expect_equal(dim(draws)[3], 2)
})
