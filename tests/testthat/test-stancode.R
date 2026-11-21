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
tdf_data_raw <- data.frame(
  Source = rep(c("Beaver", "Deer"), each = 3),
  d13C = c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3),
  d15N = c(2.9, 3.0, 3.1, 3.0, 3.1, 3.2)
)

test_that("make_stancode returns a bsimms_stancode object", {
  code <- make_stancode(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )
  expect_s3_class(code, "bsimms_stancode")
  expect_type(code, "character")
})

test_that("print.bsimms_stancode prints plain text without quotes", {
  code <- make_stancode(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = "d13C"
  )
  out <- capture.output(print(code))
  expect_false(any(grepl('^\\[1\\] "', out)))
  expect_true(any(grepl("^functions \\{", out)))
})

test_that("ilr_global floors p_global away from the simplex boundary", {
  code <- make_stancode(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )
  expect_match(code, "fmax(gmean(p_global, d), 1e-12)", fixed = TRUE)
  expect_match(code, "fmax(p_global[d + 1], 1e-12)", fixed = TRUE)
})

test_that("fixed-effect covariates add a beta parameter and per-coef priors", {
  code <- make_stancode(
    formula = ~Sex,
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )
  expect_match(code, "matrix\\[P, D\\] beta;", fixed = FALSE)
  expect_match(code, "beta\\[1\\] ~ normal\\(0, 1\\);  // SexM")
})

test_that("an intercept-only formula has no beta parameter and no P/X in data", {
  code <- make_stancode(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )
  expect_no_match(code, "matrix\\[P, D\\] beta;")
  expect_no_match(code, "int<lower=1> P;")
  expect_no_match(code, "matrix\\[N, P\\] X;")
})

test_that("a fixed-effect formula declares P/X in the data block", {
  code <- make_stancode(
    formula = ~Sex,
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )
  expect_match(code, "int<lower=1> P;")
  expect_match(code, "matrix\\[N, P\\] X;")
})

test_that("group-level terms declare sd_re_*/z_re_* and are omitted without them", {
  code_re <- make_stancode(
    formula = ~ 1 + (1 | Region),
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )
  expect_match(code_re, "sd_re_Region")
  expect_match(code_re, "z_re_Region")

  code_no_re <- make_stancode(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )
  expect_no_match(code_no_re, "sd_re_")
})

test_that("raw source data declares source_raw/source_idx and a raw-source likelihood", {
  code <- make_stancode(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = "d13C"
  )
  expect_match(code, "matrix\\[N_source_raw, J\\] source_raw;")
  expect_match(
    code,
    "source_raw\\[n\\] ~ normal\\(source_mean\\[source_idx\\[n\\]\\], source_sd\\[source_idx\\[n\\]\\]\\);"
  )
})

test_that("raw TDF data declares tdf_raw/tdf_idx, tdf_mean/tdf_sd parameters, and a raw-TDF likelihood", {
  code <- make_stancode(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_summary,
    tdf_data = tdf_data_raw,
    isotope_names = "d13C",
    source_means_sds = TRUE,
    tdf_means_sds = FALSE
  )
  expect_match(code, "int<lower=1> N_tdf_raw;")
  expect_match(
    code,
    "array\\[N_tdf_raw\\] int<lower=1, upper=K> tdf_idx;"
  )
  expect_match(code, "matrix\\[N_tdf_raw, J\\] tdf_raw;")
  expect_match(
    code,
    "matrix\\[K, J\\] tdf_mean;  // estimated TDF isotope means \\(raw TDF data\\)"
  )
  expect_match(
    code,
    "matrix<lower=0>\\[K, J\\] tdf_sd;  // estimated TDF isotope SDs \\(raw TDF data\\)"
  )
  expect_match(code, "tdf_mean\\[1, 1\\] ~ ")
  expect_match(code, "tdf_sd\\[1, 1\\] ~ ")
  expect_match(
    code,
    "tdf_raw\\[n\\] ~ normal\\(tdf_mean\\[tdf_idx\\[n\\]\\], tdf_sd\\[tdf_idx\\[n\\]\\]\\);"
  )
})

test_that("summarised source data aliases source_mean_data/source_sd_data instead", {
  code <- make_stancode(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_summary,
    tdf_data = tdf_data_summary,
    isotope_names = "d13C",
    source_means_sds = TRUE
  )
  expect_match(code, "matrix\\[K, J\\] source_mean_data;")
  expect_match(code, "matrix\\[K, J\\] source_mean = source_mean_data;")
  expect_no_match(code, "source_raw")
})

test_that("residual_only declares sigma and uses a per-isotope normal likelihood", {
  code <- make_stancode(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_summary,
    tdf_data = tdf_data_summary,
    isotope_names = "d13C",
    source_means_sds = TRUE,
    error_structure = "residual_only"
  )
  expect_match(code, "vector<lower=0>\\[J\\] sigma;")
  expect_match(
    code,
    "to_vector\\(y\\[i\\]\\) ~ normal\\(to_vector\\(mu\\[i\\]\\), sigma\\);"
  )
  expect_no_match(code, "resid_prop")
})

test_that("process_residual declares resid_prop and process_only does not", {
  code_pr <- make_stancode(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_summary,
    tdf_data = tdf_data_summary,
    isotope_names = "d13C",
    source_means_sds = TRUE,
    error_structure = "process_residual"
  )
  expect_match(code_pr, "resid_prop")
  expect_match(code_pr, "vector<lower=0, upper=20>\\[J\\] resid_prop;")

  code_po <- make_stancode(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_summary,
    tdf_data = tdf_data_summary,
    isotope_names = "d13C",
    source_means_sds = TRUE,
    error_structure = "process_only"
  )
  expect_no_match(code_po, "resid_prop")
  expect_match(
    code_po,
    "y\\[i, j\\] ~ normal\\(mu\\[i, j\\], sqrt\\(proc_var\\[i, j\\]\\)\\);"
  )
})

test_that("cross-tracer covariance appears only with raw source data and 2+ isotopes", {
  code_cov <- make_stancode(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )
  expect_match(code_cov, "Lcorr_source")
  expect_match(code_cov, "multi_normal_cholesky")

  code_single <- make_stancode(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = "d13C"
  )
  expect_no_match(code_single, "Lcorr_source")
})

test_that("residual-error correlation appears only under residual_only with 2+ isotopes", {
  code <- make_stancode(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_summary,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N"),
    source_means_sds = TRUE,
    error_structure = "residual_only"
  )
  expect_match(code, "Lcorr_resid")
  expect_match(
    code,
    "multi_normal_cholesky_lpdf\\(to_vector\\(y\\[i\\]\\) \\| to_vector\\(mu\\[i\\]\\), L_Sigma_resid\\)"
  )
})

test_that("conc_dep = TRUE declares the concentration matrix and weights proportions by it", {
  d <- source_data_raw
  d$d13C_conc <- c(0.3, 0.32, 0.28, 0.35, 0.33, 0.37)
  d$d15N_conc <- c(0.05, 0.05, 0.06, 0.04, 0.05, 0.04)
  code <- make_stancode(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = d,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N"),
    conc_dep = TRUE
  )
  expect_match(code, "matrix<lower=0, upper=1>\\[K, J\\] conc;")
  expect_match(code, "concentration-weighted")
})

test_that("center/scale standardise only continuous columns in transformed data", {
  mixture_data_x <- cbind(mixture_data, Age = c(20, 30, 40, 50, 60, 70))

  code_both <- make_stancode(
    formula = ~ Age + Sex,
    mixture_data = mixture_data_x,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )
  expect_match(code_both, "means_X\\[1\\] = mean\\(X\\[, 1\\]\\);  // Age")
  expect_match(code_both, "sds_X\\[1\\] = sd\\(X\\[, 1\\]\\);  // Age")
  expect_match(code_both, "means_X\\[2\\] = 0;  // SexM")
  expect_match(code_both, "sds_X\\[2\\] = 1;  // SexM")
  expect_match(
    code_both,
    "Xc\\[, p\\] = \\(X\\[, p\\] - means_X\\[p\\]\\) / sds_X\\[p\\];"
  )
  expect_match(code_both, "Xc * beta", fixed = TRUE)
  expect_match(code_both, "matrix\\[P, D\\] beta_raw;")
  expect_match(code_both, "vector\\[D\\] ilr_global_raw;")
  expect_match(code_both, "// covariates: centred and scaled", fixed = TRUE)

  code_center_only <- make_stancode(
    formula = ~ Age + Sex,
    mixture_data = mixture_data_x,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N"),
    scale = FALSE
  )
  expect_no_match(code_center_only, "sds_X")
  expect_no_match(code_center_only, "matrix\\[P, D\\] beta_raw;")
  expect_match(code_center_only, "vector\\[D\\] ilr_global_raw;")
  expect_match(code_center_only, "Xc\\[, p\\] = X\\[, p\\] - means_X\\[p\\];")
  expect_match(code_center_only, "// covariates: centred\n", fixed = TRUE)

  code_scale_only <- make_stancode(
    formula = ~ Age + Sex,
    mixture_data = mixture_data_x,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N"),
    center = FALSE
  )
  expect_no_match(code_scale_only, "means_X")
  expect_no_match(code_scale_only, "vector\\[D\\] ilr_global_raw;")
  expect_match(code_scale_only, "matrix\\[P, D\\] beta_raw;")
  expect_match(code_scale_only, "Xc\\[, p\\] = X\\[, p\\] / sds_X\\[p\\];")
  expect_match(code_scale_only, "// covariates: scaled\n", fixed = TRUE)

  code_neither <- make_stancode(
    formula = ~ Age + Sex,
    mixture_data = mixture_data_x,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N"),
    center = FALSE,
    scale = FALSE
  )
  expect_no_match(code_neither, "means_X|sds_X|Xc|beta_raw|ilr_global_raw")
  expect_no_match(code_neither, "// covariates:")
  expect_match(code_neither, "X * beta", fixed = TRUE)
})

test_that("center/scale generate no standardisation code for an all-factor formula", {
  code <- make_stancode(
    formula = ~Sex,
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )
  expect_no_match(code, "means_X|sds_X|Xc|beta_raw|ilr_global_raw")
  expect_no_match(code, "// covariates:")
  expect_match(code, "X * beta", fixed = TRUE)
})

test_that("indent prefixes lines with spaces, and leaves an empty vector empty", {
  expect_equal(indent(c("a", "b")), c("  a", "  b"))
  expect_equal(indent(c("a", "b"), n = 4), c("    a", "    b"))
  expect_equal(indent(character(0)), character(0))
})

test_that("standardisation_label reflects has_ilr_global_raw/has_beta_raw", {
  spec <- function(ilr_global_raw, beta_raw) {
    list(has_ilr_global_raw = ilr_global_raw, has_beta_raw = beta_raw)
  }
  expect_equal(standardisation_label(spec(TRUE, TRUE)), "centred and scaled")
  expect_equal(standardisation_label(spec(TRUE, FALSE)), "centred")
  expect_equal(standardisation_label(spec(FALSE, TRUE)), "scaled")
  expect_null(standardisation_label(spec(FALSE, FALSE)))
})

test_that("a user prior override replaces the default in the generated code", {
  code_default <- make_stancode(
    formula = ~Sex,
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )
  expect_match(code_default, "beta\\[1\\] ~ normal\\(0, 1\\);")

  code_override <- make_stancode(
    formula = ~Sex,
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N"),
    prior = bsimms_prior(prior = "normal(0, 5)", class = "b")
  )
  expect_match(code_override, "beta\\[1\\] ~ normal\\(0, 5\\);")
})

test_that("generated Stan code is syntactically valid across error structures/data modes", {
  skip_if_not_installed("cmdstanr")
  check_syntax <- function(code) {
    f <- tempfile(fileext = ".stan")
    on.exit(unlink(f))
    writeLines(as.character(code), f)
    mod <- cmdstanr::cmdstan_model(f, compile = FALSE)
    mod$check_syntax(quiet = TRUE)
  }

  expect_true(check_syntax(make_stancode(
    formula = ~ Sex + (1 | Region),
    mixture_data = mixture_data,
    source_data = source_data_raw,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N")
  )))
  expect_true(check_syntax(make_stancode(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_summary,
    tdf_data = tdf_data_summary,
    isotope_names = c("d13C", "d15N"),
    source_means_sds = TRUE,
    error_structure = "residual_only"
  )))
  expect_true(check_syntax(make_stancode(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_summary,
    tdf_data = tdf_data_summary,
    isotope_names = "d13C",
    source_means_sds = TRUE,
    error_structure = "process_only"
  )))
  expect_true(check_syntax(make_stancode(
    formula = ~1,
    mixture_data = mixture_data,
    source_data = source_data_summary,
    tdf_data = tdf_data_raw,
    isotope_names = c("d13C", "d15N"),
    source_means_sds = TRUE,
    tdf_means_sds = FALSE
  )))

  mixture_data_x <- cbind(mixture_data, Age = c(20, 30, 40, 50, 60, 70))
  for (ctr in c(TRUE, FALSE)) {
    for (scl in c(TRUE, FALSE)) {
      expect_true(check_syntax(make_stancode(
        formula = ~ Age + Sex,
        mixture_data = mixture_data_x,
        source_data = source_data_raw,
        tdf_data = tdf_data_summary,
        isotope_names = c("d13C", "d15N"),
        center = ctr,
        scale = scl
      )))
    }
  }
})
