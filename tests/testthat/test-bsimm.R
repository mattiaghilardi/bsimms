mixture_data <- data.frame(
  d13C = c(-20, -21, -19, -22, -20.5, -21.5),
  d15N = c(10, 11, 9, 12, 10.5, 11.5),
  Region = factor(rep(c("A", "B"), each = 3))
)
source_data <- data.frame(
  Source = c("Beaver", "Deer"),
  d13C_mean = c(-25, -18),
  d13C_sd = c(1, 1),
  d15N_mean = c(5, 8),
  d15N_sd = c(1, 1)
)
tdf_data <- data.frame(
  Source = c("Beaver", "Deer"),
  d13C_mean = c(1, 1.2),
  d13C_sd = c(0.2, 0.3),
  d15N_mean = c(3, 3.1),
  d15N_sd = c(0.4, 0.5)
)

test_that("resolve_backend resolves auto to an installed backend", {
  skip_if_not_installed("cmdstanr")
  expect_equal(resolve_backend("auto"), "cmdstanr")
})

test_that("resolve_backend validates an explicit backend is installed", {
  skip_if_not_installed("rstan")
  expect_equal(resolve_backend("rstan"), "rstan")
})

test_that("bsimm errors on an invalid error_structure or backend", {
  expect_snapshot(
    error = TRUE,
    bsimm(
      formula = ~1,
      mixture_data = mixture_data,
      source_data = source_data,
      tdf_data = tdf_data,
      isotope_names = c("d13C", "d15N"),
      error_structure = "banana",
      source_means_sds = TRUE
    )
  )
  expect_snapshot(
    error = TRUE,
    bsimm(
      formula = ~1,
      mixture_data = mixture_data,
      source_data = source_data,
      tdf_data = tdf_data,
      isotope_names = c("d13C", "d15N"),
      backend = "banana",
      source_means_sds = TRUE
    )
  )
})

test_that("bsimm fits via cmdstanr and returns a well-formed bsimms object", {
  skip_if_not_installed("cmdstanr")
  fit <- bsimm(
    formula = ~ 1 + (1 | Region),
    mixture_data = mixture_data,
    source_data = source_data,
    tdf_data = tdf_data,
    isotope_names = c("d13C", "d15N"),
    source_means_sds = TRUE,
    chains = 1,
    iter_warmup = 200,
    iter_sampling = 100,
    seed = 1,
    refresh = 0
  )
  expect_s3_class(fit, "bsimms_fit")
  expect_equal(fit$backend, "cmdstanr")
  expect_s3_class(fit$fit, "CmdStanMCMC")
  expect_s3_class(fit$stancode, "bsimms_stancode")
  expect_s3_class(fit$standata, "bsimms_standata")
  expect_s3_class(fit$prior, "bsimms_prior")
  expect_s3_class(fit$spec, "bsimms_spec")
  expect_equal(fit$call[[1]], as.name("bsimm"))

  draws <- fit$fit$draws("p_global", format = "matrix")
  expect_equal(colnames(draws), c("p_global[1]", "p_global[2]"))
  expect_true(all(rowSums(draws) > 0.99 & rowSums(draws) < 1.01))
})

test_that("bsimm fits via rstan when explicitly requested", {
  skip_if_not_installed("rstan")
  skip_on_os("windows") # rstan model compilation is unreliably slow/fragile in Windows CI
  fit <- suppressWarnings(bsimm(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data,
    tdf_data = tdf_data,
    isotope_names = c("d13C", "d15N"),
    backend = "rstan",
    source_means_sds = TRUE,
    chains = 1,
    iter_warmup = 200,
    iter_sampling = 100,
    seed = 1,
    refresh = 0
  ))
  expect_equal(fit$backend, "rstan")
  expect_s4_class(fit$fit, "stanfit")
})

test_that("a user prior override is reflected in the fitted object", {
  skip_if_not_installed("cmdstanr")
  fit <- bsimm(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data,
    tdf_data = tdf_data,
    isotope_names = c("d13C", "d15N"),
    source_means_sds = TRUE,
    prior = bsimms_prior(prior = "3", class = "p_global", group = "Beaver"),
    chains = 1,
    iter_warmup = 200,
    iter_sampling = 100,
    seed = 1,
    refresh = 0
  )
  expect_equal(
    fit$prior$prior[
      fit$prior$class == "p_global" & fit$prior$group == "Beaver"
    ],
    "3"
  )
  expect_equal(fit$standata$alpha_dirichlet, c(3, 1))
})
