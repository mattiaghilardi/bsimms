# error_structure defaults to process_residual and is validated

    Code
      build_bsimms_spec(formula = ~1, mixture_data = mixture_data, source_data = source_data,
        tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), error_structure = "banana")
    Condition
      Error in `build_bsimms_spec()`:
      ! `error_structure` must be one of "process_residual", "process_only", or "residual_only", not "banana".

# a single mixture data point requires process_only error

    Code
      build_bsimms_spec(formula = ~1, mixture_data = one_mix, source_data = source_data,
        tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), error_structure = "process_residual")
    Condition
      Error:
      ! With a single mixture data point, `error_structure` must be "process_only".
      i Residual (observation-level) error cannot be estimated without replication; only process error, propagated from the source/TDF data, is identifiable here.

---

    Code
      build_bsimms_spec(formula = ~1, mixture_data = one_mix, source_data = source_data,
        tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), error_structure = "residual_only")
    Condition
      Error:
      ! With a single mixture data point, `error_structure` must be "process_only".
      i Residual (observation-level) error cannot be estimated without replication; only process error, propagated from the source/TDF data, is identifiable here.

# one mixture sample per level of a single fixed factor requires process_only

    Code
      build_bsimms_spec(formula = ~Individual, mixture_data = ind_mix, source_data = source_data,
        tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), error_structure = "process_residual")
    Condition
      Error:
      ! With one mixture data point per level of Individual, `error_structure` must be "process_only".
      i Residual (observation-level) error cannot be separated from the effect of Individual when every level has exactly one mixture sample.

# one mixture sample per level of a random-effect group requires process_only

    Code
      build_bsimms_spec(formula = ~ 1 + (1 | Individual), mixture_data = ind_mix,
      source_data = source_data, tdf_data = tdf_data, isotope_names = c("d13C",
        "d15N"), error_structure = "process_residual")
    Condition
      Error:
      ! With one mixture data point per level of Individual, `error_structure` must be "process_only".
      i Residual (observation-level) error cannot be separated from the group-level variance of Individual when every level has exactly one mixture sample.

# build_bsimms_spec errors on invalid isotope_names

    Code
      build_bsimms_spec(formula = ~1, mixture_data = mixture_data, source_data = source_data,
        tdf_data = tdf_data, isotope_names = character(0))
    Condition
      Error:
      ! `isotope_names` must be a non-empty character vector of mixture/source/TDF column names.

# build_bsimms_spec errors with fewer than 2 sources

    Code
      build_bsimms_spec(formula = ~1, mixture_data = mixture_data, source_data = one_source,
        tdf_data = one_tdf, isotope_names = c("d13C", "d15N"))
    Condition
      Error:
      ! Need at least 2 sources.

# build_bsimms_spec errors when the formula's design matrix drops rows

    Code
      build_bsimms_spec(formula = ~Region, mixture_data = d, source_data = source_data,
        tdf_data = tdf_data, isotope_names = c("d13C", "d15N"))
    Condition
      Error:
      ! Design matrix built from `formula` has 5 rows, but the mixture isotope data has 6.
      i This is usually caused by missing values in the fixed-effect covariates: rows with `NA` are silently dropped when building the design matrix.
      i Remove or impute the missing covariate values in `mixture_data` before fitting.

