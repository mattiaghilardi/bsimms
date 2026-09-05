#' Simulate data for a `bsimms` model
#'
#' Generates `mixture_data`, `source_data`, and `tdf_data` for an arbitrary
#' `bsimms` model configuration (any number of sources/isotopes, raw or
#' summarised source/TDF data, any error structure, with or without
#' concentration dependence, and any `lme4`-style fixed/random-effects
#' `formula` of mixture-level covariates), together with the known
#' generating ("true") parameter values. Useful for testing `bsimm()`
#' across configurations, and for checking parameter recovery.
#'
#' Source proportions and mixture isotope values are simulated forward from
#' known truth: a baseline `p_global`, formula-driven fixed/random-effect
#' deviations from it in ILR space, source and TDF isotope means/SDs, and
#' observation-level noise appropriate to `error_structure`. Random-effect
#' terms in `formula` are currently limited to random intercepts (`(1 |
#' Group)`, including crossed/nested grouping factors); random slopes are
#' not yet supported.
#'
#' Every covariate referenced in `formula` must be either a fixed-effect
#' factor named in `n_levels`, a random-effect grouping factor named in
#' `n_groups`, or otherwise is generated as an independent standard normal
#' (continuous) covariate. Source means are generated to be well-separated
#' and non-collinear in isotope space (sampled without replacement from a
#' pool of `2 * n_sources` equally spaced candidate values per isotope,
#' redrawn if the resulting configuration is collinear), rather than
#' sampled independently from a probability distribution, which cannot
#' guarantee isotopically distinct sources.
#'
#' A warning is issued when `n_sources > n_isotopes + 1`: the mixing system
#' is then underdetermined in the classical sense (more sources than
#' isotopes can resolve without help from the prior), which `bsimms`
#' supports but which relies more heavily on the prior than the data.
#'
#' @param formula A one-sided `lme4`-style formula, as in [bsimm()]. Random
#'   slopes (e.g. `(x | Group)`) are not supported.
#' @param n_mixture_obs Number of mixture observations to simulate.
#' @param n_sources Number of sources to simulate. Ignored if `source_names`
#'   is supplied.
#' @param source_names Optional character vector of source names. If
#'   `NULL` (default), `n_sources` sources named `"source1"`, `"source2"`,
#'   ... are used.
#' @param n_isotopes Number of isotopes to simulate. Ignored if
#'   `isotope_names` is supplied.
#' @param isotope_names Optional character vector of isotope names. If
#'   `NULL` (default), `n_isotopes` isotopes named `"isotope1"`,
#'   `"isotope2"`, ... are used.
#' @param n_levels Named list, one entry per fixed-effect factor in
#'   `formula`, giving that factor's number of levels (an integer >= 2).
#'   Any variable in `formula` not listed here is generated as a continuous
#'   covariate instead.
#' @param n_groups Named list, one entry per random-effect grouping factor
#'   in `formula` (i.e. every variable referenced to the right of `|`),
#'   giving that factor's number of groups (an integer >= 2). Required for
#'   every grouping factor in `formula`.
#' @param balanced Logical; split `n_mixture_obs` as evenly as possible
#'   across each factor's/grouping factor's levels (`TRUE`, default; sizes
#'   differ by at most one observation when `n_mixture_obs` is not a
#'   multiple of the number of levels), or assign them at random (`FALSE`)?
#'   Either way every level/group gets at least one observation.
#' @param source_means_sds,tdf_means_sds Logical; generate `source_data`/
#'   `tdf_data` as means/SDs (`TRUE`) or raw replicate samples (`FALSE`)?
#'   Default `FALSE` for `source_means_sds`, `TRUE` for `tdf_means_sds`,
#'   matching [bsimm()]'s own defaults.
#' @param conc_dep Logical; also generate a `<isotope>_conc` column (in
#'   `(0, 1]`, summing to at most 1 across isotopes for a given source) in
#'   `source_data` for every isotope, enabling concentration dependence
#'   (default `FALSE`)?
#' @inheritParams bsimm
#' @param n_source_obs Number of raw replicate observations to simulate per
#'   source: either a single integer, recycled to every source (default
#'   10), or a numeric vector of length `n_sources` giving each source its
#'   own count (allowing unbalanced source data). Ignored if
#'   `source_means_sds = TRUE`.
#' @param n_tdf_obs As `n_source_obs`, for `tdf_data`. Ignored if
#'   `tdf_means_sds = TRUE`.
#' @param p_global Optional numeric vector of length `n_sources`, the true
#'   baseline source proportions (must sum to 1). If `NULL` (default),
#'   randomly generated.
#' @param sigma Optional numeric vector of length `n_isotopes`, the true
#'   residual SD for each isotope. Used only when `error_structure =
#'   "residual_only"` (ignored otherwise). If `NULL` (default), randomly
#'   generated.
#' @param resid_prop Optional numeric vector of length `n_isotopes`, the
#'   true residual-error factor scaling process variance for each isotope.
#'   Used only when `error_structure = "process_residual"` (ignored
#'   otherwise). If `NULL` (default), randomly generated.
#' @param seed Optional integer seed. The global RNG state is saved and
#'   restored, so calling this function does not affect the caller's own
#'   random draws.
#' @return A list with elements `formula`, `mixture_data`, `source_data`,
#'   `tdf_data`, `isotope_names`, `source_means_sds`, `tdf_means_sds`,
#'   `conc_dep`, `error_structure`, `source_col` (in the same order as
#'   [bsimm()]'s arguments, ready to pass to it), and `truth`, the true
#'   generating values that [bsimm()] would later estimate from the data:
#'   `p_global`, `fixed` (matrix of fixed-effect ILR-space coefficients, if
#'   `formula` has fixed-effect terms), `random` (list of group-level
#'   SDs/realised effects, one element per random-effect term, if any),
#'   `source_mean`/`source_sd` (if `source_means_sds = FALSE`),
#'   `tdf_mean`/`tdf_sd` (if `tdf_means_sds = FALSE`), and `sigma` or
#'   `resid_prop` (whichever applies to `error_structure`).
#' @inherit make_stancode references
#' @export
#' @examples
#' sim <- simulate_bsimms_data(
#'   ~ Sex + (1 | Region),
#'   n_mixture_obs = 60,
#'   n_levels = list(Sex = 2),
#'   n_groups = list(Region = 3)
#' )
#' str(sim, max.level = 1)
simulate_bsimms_data <- function(
  formula = ~1,
  n_mixture_obs,
  n_sources = 3,
  source_names = NULL,
  n_isotopes = 2,
  isotope_names = NULL,
  n_levels = list(),
  n_groups = list(),
  balanced = TRUE,
  source_means_sds = FALSE,
  n_source_obs = 10,
  tdf_means_sds = TRUE,
  n_tdf_obs = 10,
  conc_dep = FALSE,
  error_structure = c("process_residual", "process_only", "residual_only"),
  p_global = NULL,
  sigma = NULL,
  resid_prop = NULL,
  source_col = "Source",
  seed = NULL
) {
  error_structure <- rlang::arg_match(error_structure)
  if (!inherits(formula, "formula") || length(formula) != 2) {
    cli::cli_abort(
      paste0(
        "{.arg formula} must be a one-sided formula, ",
        "e.g. {.code ~ Sex + (1 | Region)}."
      ),
      call = NULL
    )
  }
  if (!rlang::is_scalar_integerish(n_mixture_obs) || n_mixture_obs < 1) {
    cli::cli_abort(
      "{.arg n_mixture_obs} must be a positive integer.",
      call = NULL
    )
  }

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (!is.null(old_seed)) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    })
    set.seed(seed)
  }

  source_names <- if (is.null(source_names)) {
    sim_default_names(n_sources, "source")
  } else {
    source_names
  }
  K <- length(source_names)
  isotope_names <- if (is.null(isotope_names)) {
    sim_default_names(n_isotopes, "isotope")
  } else {
    isotope_names
  }
  J <- length(isotope_names)
  D <- K - 1L
  n_source_obs <- sim_resolve_n_obs(n_source_obs, K, "n_source_obs")
  n_tdf_obs <- sim_resolve_n_obs(n_tdf_obs, K, "n_tdf_obs")

  if (K > J + 1) {
    cli::cli_warn(
      c(
        paste0(
          "{K} sources with only {J} isotope{?s} is an underdetermined ",
          "mixing system ({.code n_sources > n_isotopes + 1})."
        ),
        "i" = paste0(
          "Source proportions will rely more heavily on the prior than ",
          "the data."
        )
      )
    )
  }

  # Parse `formula` just enough to know which base variables need to be
  # generated, without needing data yet -- `parse_bsimms_formula()` builds
  # the actual design matrices once those columns exist (see below).
  expanded <- reformulas::expandDoubleVerts(formula)
  bars <- reformulas::findbars(expanded)
  fixed_formula <- reformulas::nobars(expanded)
  fixed_vars <- all.vars(fixed_formula)

  group_vars <- character(0)
  for (b in bars) {
    b_str <- rlang::expr_deparse(b)
    parts <- strsplit(b_str, "\\|")[[1]]
    term_str <- trimws(parts[1])
    group_str <- trimws(parts[2])
    if (term_str != "1") {
      cli::cli_abort(
        c(
          "{.arg formula} contains a random slope ({.code {b_str}}).",
          "i" = paste0(
            "{.fn simulate_bsimms_data} only supports random intercepts, ",
            "e.g. {.code (1 | Group)}."
          )
        ),
        call = NULL
      )
    }
    group_vars <- c(group_vars, all.vars(str2lang(group_str)))
  }
  group_vars <- unique(group_vars)

  sim_check_condition_names(n_levels, "n_levels", fixed_vars)
  sim_check_condition_names(n_groups, "n_groups", group_vars)
  missing_groups <- setdiff(group_vars, names(n_groups))
  if (length(missing_groups) > 0) {
    cli::cli_abort(
      paste0(
        "{.arg n_groups} must include an entry for every grouping ",
        "factor in {.arg formula}: missing {.val {missing_groups}}."
      ),
      call = NULL
    )
  }

  mixture_data <- data.frame(row.names = seq_len(n_mixture_obs))
  for (v in fixed_vars) {
    mixture_data[[v]] <- if (v %in% names(n_levels)) {
      sim_factor_column(
        n_mixture_obs,
        n_levels[[v]],
        balanced,
        sim_level_names(n_levels[[v]])
      )
    } else {
      stats::rnorm(n_mixture_obs)
    }
  }
  for (v in setdiff(group_vars, names(mixture_data))) {
    mixture_data[[v]] <- sim_factor_column(
      n_mixture_obs,
      n_groups[[v]],
      balanced,
      sim_level_names(n_groups[[v]])
    )
  }

  parsed <- parse_bsimms_formula(formula, data = mixture_data)
  P <- length(parsed$fixed_names) -
    as.integer("(Intercept)" %in% parsed$fixed_names)
  X <- parsed$X[, parsed$fixed_names != "(Intercept)", drop = FALSE]

  # Truth: baseline + fixed/random-effect deviations, in ILR space
  if (is.null(p_global)) {
    g <- stats::rgamma(K, shape = 2)
    p_global <- g / sum(g)
  } else if (
    length(p_global) != K ||
      !isTRUE(all.equal(sum(p_global), 1)) ||
      any(p_global <= 0)
  ) {
    cli::cli_abort(
      paste0(
        "{.arg p_global} must be a vector of {K} strictly positive ",
        "proportions summing to 1."
      ),
      call = NULL
    )
  }
  eta <- matrix(
    rep(ilr(p_global), each = n_mixture_obs),
    nrow = n_mixture_obs,
    ncol = D
  )

  truth <- list(p_global = p_global)

  beta <- NULL
  if (P > 0) {
    beta <- matrix(
      stats::rnorm(P * D, 0, 0.5),
      nrow = P,
      ncol = D,
      dimnames = list(colnames(X), NULL)
    )
    eta <- eta + X %*% beta
    truth$fixed <- beta
  }

  truth$random <- list()
  for (re in parsed$re_terms) {
    n_grp <- length(re$group_levels)
    sd_re <- abs(stats::rnorm(D, mean = 0.3, sd = 0.1))
    b_re <- matrix(
      stats::rnorm(D * n_grp, 0, rep(sd_re, times = n_grp)),
      nrow = D,
      ncol = n_grp
    )
    dimnames(b_re) <- list(NULL, re$group_levels)
    eta <- eta + t(b_re[, re$group_idx, drop = FALSE])
    truth$random[[re$label]] <- list(
      group = re$group,
      sd = sd_re,
      effects = b_re
    )
  }
  if (length(truth$random) == 0) {
    truth$random <- NULL
  }

  p_true <- ilr_inv(eta)

  # Truth: source and TDF isotope means/SDs
  source_mean <- sim_separated_source_means(K, J)
  dimnames(source_mean) <- list(source_names, isotope_names)
  source_sd <- matrix(
    stats::runif(K * J, 0.5, 1.5),
    K,
    J,
    dimnames = dimnames(source_mean)
  )

  tdf_mean <- matrix(
    2 + stats::rnorm(K * J, 0, 0.2),
    K,
    J,
    dimnames = dimnames(source_mean)
  )
  tdf_sd <- matrix(
    stats::runif(K * J, 0.3, 0.7),
    K,
    J,
    dimnames = dimnames(source_mean)
  )

  if (!source_means_sds) {
    truth$source_mean <- source_mean
  }
  if (!source_means_sds) {
    truth$source_sd <- source_sd
  }
  if (!tdf_means_sds) {
    truth$tdf_mean <- tdf_mean
  }
  if (!tdf_means_sds) {
    truth$tdf_sd <- tdf_sd
  }

  conc <- NULL
  if (conc_dep) {
    conc <- matrix(0, K, J, dimnames = dimnames(source_mean))
    for (k in seq_len(K)) {
      w <- stats::rgamma(J, shape = 2)
      conc[k, ] <- w / sum(w) * stats::runif(1, 0.3, 0.9)
    }
  }

  needs_proc <- error_structure != "residual_only"
  needs_resid_prop <- error_structure == "process_residual"
  if (needs_proc) {
    mu_var <- sim_mu_var(p_true, source_mean, source_sd, tdf_mean, tdf_sd, conc)
    mu <- mu_var$mu
    proc_var <- mu_var$proc_var
  } else {
    mu <- sim_mu_var(
      p_true,
      source_mean,
      source_sd,
      tdf_mean,
      tdf_sd,
      conc,
      need_var = FALSE
    )$mu
  }

  if (error_structure == "residual_only") {
    if (is.null(sigma)) {
      sigma <- stats::runif(J, 0.3, 0.6)
    } else if (length(sigma) != J) {
      cli::cli_abort(
        "{.arg sigma} must be a numeric vector of length {J}.",
        call = NULL
      )
    }
    truth$sigma <- stats::setNames(sigma, isotope_names)
    y_sd <- matrix(rep(sigma, each = n_mixture_obs), n_mixture_obs, J)
  } else if (error_structure == "process_residual") {
    if (is.null(resid_prop)) {
      resid_prop <- stats::runif(J, 0.5, 3)
    } else if (length(resid_prop) != J) {
      cli::cli_abort(
        "{.arg resid_prop} must be a numeric vector of length {J}.",
        call = NULL
      )
    }
    truth$resid_prop <- stats::setNames(resid_prop, isotope_names)
    y_sd <- sqrt(sweep(proc_var, 2, resid_prop, `*`))
  } else {
    y_sd <- sqrt(proc_var)
  }
  y <- mu + matrix(stats::rnorm(n_mixture_obs * J, 0, y_sd), n_mixture_obs, J)
  colnames(y) <- isotope_names
  mixture_data <- cbind(mixture_data, as.data.frame(y))

  source_data <- sim_isotope_data(
    source_names,
    isotope_names,
    source_mean,
    source_sd,
    source_means_sds,
    n_source_obs,
    source_col,
    conc
  )
  tdf_data <- sim_isotope_data(
    source_names,
    isotope_names,
    tdf_mean,
    tdf_sd,
    tdf_means_sds,
    n_tdf_obs,
    source_col
  )

  list(
    formula = formula,
    mixture_data = mixture_data,
    source_data = source_data,
    tdf_data = tdf_data,
    isotope_names = isotope_names,
    source_means_sds = source_means_sds,
    tdf_means_sds = tdf_means_sds,
    conc_dep = conc_dep,
    error_structure = error_structure,
    source_col = source_col,
    truth = truth
  )
}

#' Default `"source1"`/`"isotope1"`-style names, zero-padded so lexicographic
#' sorting (as done downstream by `build_bsimms_spec()`) matches numeric
#' order regardless of `n`.
#'
#' @param n Integer, number of names.
#' @param prefix Character prefix.
#' @return Character vector of length `n`.
#' @noRd
sim_default_names <- function(n, prefix) {
  sprintf(paste0(prefix, "%0", nchar(n), "d"), seq_len(n))
}

#' Default factor level names: `"A"`, `"B"`, ..., falling back to `"L1"`,
#' `"L2"`, ... beyond 26 levels.
#'
#' @param n Integer, number of levels.
#' @return Character vector of length `n`.
#' @noRd
sim_level_names <- function(n) {
  if (n <= 26) LETTERS[seq_len(n)] else paste0("L", seq_len(n))
}

#' Validate that every name in a `n_levels`/`n_groups`-style named list is
#' an actual formula variable, and that every count is a valid (>= 2)
#' number of levels.
#'
#' @param x Named list.
#' @param arg_name Character, for error messages.
#' @param available Character vector of valid names.
#' @return Invisible `NULL`.
#' @noRd
sim_check_condition_names <- function(x, arg_name, available) {
  unknown <- setdiff(names(x), available)
  if (length(unknown) > 0) {
    cli::cli_abort(
      paste0(
        "{.arg {arg_name}} names variable{?s} not in {.arg formula}: ",
        "{.val {unknown}}."
      ),
      call = NULL
    )
  }
  bad <- vapply(
    x,
    function(n) !rlang::is_scalar_integerish(n) || n < 2,
    logical(1)
  )
  if (any(bad)) {
    cli::cli_abort(
      "{.arg {arg_name}} entries must be single integers >= 2.",
      call = NULL
    )
  }
  invisible(NULL)
}

#' Resolve `n_source_obs`/`n_tdf_obs` into a length-`K` vector of positive
#' integer replicate counts, one per source: a scalar is recycled to every
#' source, a length-`K` vector is used as-is (allowing unbalanced source
#' data), and anything else errors.
#'
#' @param n_obs A single integer or a numeric vector of length `K`.
#' @param K Positive integer, number of sources.
#' @param arg_name Character, for error messages.
#' @return Integer vector of length `K`.
#' @noRd
sim_resolve_n_obs <- function(n_obs, K, arg_name) {
  valid_length <- length(n_obs) %in% c(1, K)
  valid_values <- is.numeric(n_obs) &&
    all(n_obs == round(n_obs)) &&
    all(n_obs >= 1)
  if (!valid_length || !valid_values) {
    cli::cli_abort(
      paste0(
        "{.arg {arg_name}} must be a single positive integer or a ",
        "positive integer vector of length {K} (n_sources)."
      ),
      call = NULL
    )
  }
  if (length(n_obs) == 1) rep(as.integer(n_obs), K) else as.integer(n_obs)
}

#' Split `n_obs` observations across `n_lvl` levels, either as evenly as
#' possible (`balanced = TRUE`) or at random (`FALSE`), guaranteeing every
#' level gets at least one observation either way -- for the random case,
#' via a direct construction (reserve one observation per level, then
#' `rmultinom()` the remainder) rather than rejection sampling, which would
#' never terminate if `n_lvl > n_obs`.
#'
#' @param n_obs,n_lvl Positive integers, `n_lvl <= n_obs`.
#' @param balanced Logical.
#' @return Integer vector of length `n_lvl`, summing to `n_obs`.
#' @noRd
sim_level_counts <- function(n_obs, n_lvl, balanced) {
  if (n_lvl > n_obs) {
    cli::cli_abort(
      paste0(
        "Cannot split {n_obs} observations across {n_lvl} levels: ",
        "each level needs at least one observation."
      ),
      call = NULL
    )
  }
  if (balanced) {
    base <- n_obs %/% n_lvl
    rem <- n_obs %% n_lvl
    counts <- rep(base, n_lvl)
    if (rem > 0) {
      counts[seq_len(rem)] <- counts[seq_len(rem)] + 1L
    }
    counts
  } else {
    extra <- as.vector(stats::rmultinom(
      1,
      n_obs - n_lvl,
      prob = rep(1 / n_lvl, n_lvl)
    ))
    rep(1L, n_lvl) + extra
  }
}

#' Generate a shuffled factor column with `n_lvl` levels, `n_obs` rows.
#'
#' @param n_obs,n_lvl Positive integers.
#' @param balanced Logical, see `sim_level_counts()`.
#' @param level_names Character vector of length `n_lvl`.
#' @return A factor of length `n_obs`.
#' @noRd
sim_factor_column <- function(n_obs, n_lvl, balanced, level_names) {
  counts <- sim_level_counts(n_obs, n_lvl, balanced)
  values <- sample(rep(seq_len(n_lvl), times = counts))
  factor(level_names[values], levels = level_names)
}

#' Generate a `K x J` matrix of well-separated, non-collinear source means:
#' for each isotope, `K` values are sampled without replacement from a pool
#' of `2 * K` equally spaced candidates (twice as many as needed, so that
#' isotopes are not forced to reuse the same source-to-value pairing), and
#' the whole matrix is redrawn if it turns out collinear (checked exactly,
#' via `cor() == 1`, for two isotopes; via the smallest of the
#' `min(K - 1, J)` non-trivial singular values of the centred matrix, for
#' more than two). Collinearity is only possible, and only checked, when
#' there are at least 3 sources and 2 isotopes.
#'
#' @param K,J Integers, number of sources and isotopes.
#' @param sd_ref Reference SD used to scale grid spacing (default 1);
#'   individual sources' actual SDs vary modestly around this value, so
#'   spacing several times larger than `sd_ref` keeps sources separated by
#'   construction rather than by chance.
#' @param spacing_factor Grid spacing, as a multiple of `sd_ref` (default
#'   5).
#' @param max_tries Maximum number of redraws before giving up and
#'   returning the last (possibly still collinear) draw (default 100);
#'   given the size of the candidate pool this is expected to never be
#'   reached in practice.
#' @return A `K x J` numeric matrix.
#' @noRd
sim_separated_source_means <- function(
  K,
  J,
  sd_ref = 1,
  spacing_factor = 5,
  max_tries = 100
) {
  spacing <- spacing_factor * sd_ref
  pool <- seq(-(K - 0.5) * spacing, (K - 0.5) * spacing, by = spacing)

  draw_once <- function() {
    means <- matrix(0, nrow = K, ncol = J)
    for (j in seq_len(J)) {
      means[, j] <- sample(pool, K, replace = FALSE)
    }
    means
  }

  is_collinear <- function(means) {
    if (K < 3 || J < 2) {
      return(FALSE)
    }
    if (J == 2) {
      return(isTRUE(abs(stats::cor(means[, 1], means[, 2])) == 1))
    }
    rank_expected <- min(K - 1, J)
    sv <- svd(scale(means, center = TRUE, scale = FALSE))$d
    sv[rank_expected] < spacing * 1e-6
  }

  means <- draw_once()
  for (i in seq_len(max_tries)) {
    if (!is_collinear(means)) {
      break
    }
    means <- draw_once()
  }
  means
}

#' Compute the true expected mixture isotope values (`mu`) and, if needed,
#' the source/TDF process variance propagated into the mixture
#' (`proc_var`), for every mixture observation and isotope -- the R-level
#' equivalent of the `mu`/`proc_var` computation in the generated Stan
#' code's `transformed parameters` block (see
#' `stan_transformed_parameters_lines()` in `stancode.R`), including the
#' concentration-weighted form used when `conc` is supplied.
#'
#' @param p_true `n_mixture_obs x K` matrix of true source proportions.
#' @param source_mean,source_sd,tdf_mean,tdf_sd `K x J` matrices.
#' @param conc `K x J` matrix of concentrations, or `NULL` (no
#'   concentration dependence).
#' @param need_var Logical; compute `proc_var` (default `TRUE`)?
#' @return A list with `mu` (`n_mixture_obs x J`) and, if `need_var`,
#'   `proc_var` (same shape).
#' @noRd
sim_mu_var <- function(
  p_true,
  source_mean,
  source_sd,
  tdf_mean,
  tdf_sd,
  conc = NULL,
  need_var = TRUE
) {
  n <- nrow(p_true)
  J <- ncol(source_mean)
  mu <- matrix(0, n, J)
  proc_var <- if (need_var) matrix(0, n, J) else NULL
  total <- source_mean + tdf_mean
  var_k <- source_sd^2 + tdf_sd^2
  for (j in seq_len(J)) {
    if (!is.null(conc)) {
      w <- sweep(p_true, 2, conc[, j], `*`)
      w <- w / rowSums(w)
    } else {
      w <- p_true
    }
    mu[, j] <- as.numeric(w %*% total[, j])
    if (need_var) proc_var[, j] <- as.numeric(w^2 %*% var_k[, j])
  }
  list(mu = mu, proc_var = proc_var)
}

#' Build a `source_data`/`tdf_data`-shaped data frame from `K x J` truth
#' matrices: either the means/SDs directly (summarised layout) or `n_obs`
#' simulated raw replicate rows per source (raw layout), and the
#' `<isotope>_conc` columns if `conc` is supplied.
#'
#' @param source_names,isotope_names Character vectors.
#' @param mean_mat,sd_mat `K x J` matrices.
#' @param means_sds Logical; summarised (`TRUE`) or raw (`FALSE`) layout.
#' @param n_obs Integer vector of length `K`, raw replicate rows per
#'   source. Ignored if `means_sds`.
#' @param source_col Name of the source-identity column.
#' @param conc Optional `K x J` matrix of concentrations, added as
#'   `<isotope>_conc` columns (summarised layout only).
#' @return A data frame.
#' @noRd
sim_isotope_data <- function(
  source_names,
  isotope_names,
  mean_mat,
  sd_mat,
  means_sds,
  n_obs,
  source_col,
  conc = NULL
) {
  if (means_sds) {
    out <- data.frame(x = source_names)
    names(out) <- source_col
    for (j in seq_along(isotope_names)) {
      out[[paste0(isotope_names[j], "_mean")]] <- mean_mat[, j]
      out[[paste0(isotope_names[j], "_sd")]] <- sd_mat[, j]
    }
    if (!is.null(conc)) {
      for (j in seq_along(isotope_names)) {
        out[[paste0(isotope_names[j], "_conc")]] <- conc[, j]
      }
    }
    out
  } else {
    K <- length(source_names)
    out <- do.call(
      rbind,
      lapply(seq_len(K), function(k) {
        row <- data.frame(x = rep(source_names[k], n_obs[k]))
        names(row) <- source_col
        for (j in seq_along(isotope_names)) {
          row[[isotope_names[j]]] <- stats::rnorm(
            n_obs[k],
            mean_mat[k, j],
            sd_mat[k, j]
          )
        }
        if (!is.null(conc)) {
          for (j in seq_along(isotope_names)) {
            row[[paste0(isotope_names[j], "_conc")]] <- conc[k, j]
          }
        }
        row
      })
    )
    rownames(out) <- NULL
    out
  }
}
