#' Infer the sorted set of source names from `source_data`, cross-checked
#' against `tdf_data` if supplied (both must reference exactly the same
#' sources).
#'
#' @param source_data Source (diet item) data frame, raw or summarised; must
#'   have a `source_col` column.
#' @param tdf_data Trophic discrimination factor data frame, raw or
#'   summarised, with a `source_col` column; or `NULL` to skip cross-checking.
#' @param source_col Name of the source-identifier column shared by
#'   `source_data` and `tdf_data`.
#' @return Character vector of sorted unique source names.
#' @noRd
infer_source_names <- function(source_data, tdf_data, source_col = "Source") {
  if (!rlang::is_string(source_col)) {
    cli::cli_abort("{.arg source_col} must be a single string.", call = NULL)
  }
  if (!source_col %in% names(source_data)) {
    cli::cli_abort("Source data must have a {.field {source_col}} column.", call = NULL)
  }
  src_names <- sort(unique(as.character(source_data[[source_col]])))
  if (!is.null(tdf_data)) {
    if (!source_col %in% names(tdf_data)) {
      cli::cli_abort("TDF data must have a {.field {source_col}} column.", call = NULL)
    }
    tdf_names <- sort(unique(as.character(tdf_data[[source_col]])))
    if (!setequal(src_names, tdf_names)) {
      cli::cli_abort(c(
        "Source names in {.field source_data} and {.field tdf_data} do not match.",
        "x" = "{.field source_data}: {.val {src_names}}",
        "x" = "{.field tdf_data}: {.val {tdf_names}}"
      ), call = NULL)
    }
  }
  src_names
}

#' Extract and validate the mixture isotope value matrix: a numeric
#' `nrow(mixture_data) x length(isotope_names)` matrix, columns in
#' `isotope_names` order. Errors on missing columns or missing values.
#'
#' @param mixture_data Mixture data frame, with one column per isotope in
#'   `isotope_names`.
#' @param isotope_names Character vector of isotope column names to extract,
#'   in the order they should appear in the output matrix.
#' @return Numeric `nrow(mixture_data) x length(isotope_names)` matrix.
#' @noRd
prep_mixture_isotopes <- function(mixture_data, isotope_names) {
  missing_cols <- setdiff(isotope_names, names(mixture_data))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      "Mixture data is missing isotope column{?s}: {.field {missing_cols}}.",
      call = NULL
    )
  }
  y <- as.matrix(mixture_data[, isotope_names, drop = FALSE])
  storage.mode(y) <- "double"
  if (anyNA(y)) {
    cli::cli_abort("Mixture isotope data contain missing values; remove or impute them first.", call = NULL)
  }
  y
}

#' Prepare source or TDF isotope data, raw or summarised
#'
#' Internal helper shared by source and TDF data preparation. When
#' `means_sds = FALSE`, `data` is expected in long format with one row per
#' raw measurement and columns `source_col`, `isotope_names`. When
#' `means_sds = TRUE`, `data` is expected in wide format with one row per
#' source and columns `source_col`, `<isotope>_mean`, `<isotope>_sd` for
#' every isotope in `isotope_names`.
#'
#' @param data Source or TDF data frame, raw or summarised (see above).
#' @param isotope_names Character vector of isotope names.
#' @param source_names Character vector of the full set of source names (as
#'   returned by `infer_source_names()`); `data` must reference a subset of
#'   these (or, when `means_sds = TRUE`, exactly these).
#' @param means_sds Logical; is `data` summarised (means/SDs, `TRUE`) or raw
#'   replicate measurements (`FALSE`)?
#' @param source_col Name of the source-identifier column in `data`.
#' @param label Either `"source"` or `"tdf"`, used only to name the data set
#'   in error messages (e.g. `"` source_data`"`).
#' @return A list. If `means_sds = FALSE`: `mode = "raw"`, `n` (row count),
#'   `source_idx` (integer vector mapping each row to its index in
#'   `source_names`), `Y` (raw isotope value matrix). If `means_sds = TRUE`:
#'   `mode = "summary"`, `mean` and `sd` (`length(source_names) x
#'   length(isotope_names)` matrices, rows in `source_names` order).
#' @noRd
prep_iso_table <- function(data, isotope_names, source_names, means_sds,
                            source_col = "Source", label = "source") {
  label <- rlang::arg_match0(label, c("source", "tdf"))
  if (!source_col %in% names(data)) {
    cli::cli_abort("{.field {label}_data} must have a {.field {source_col}} column.", call = NULL)
  }
  data[[source_col]] <- as.character(data[[source_col]])
  unknown <- setdiff(unique(data[[source_col]]), source_names)
  if (length(unknown) > 0) {
    cli::cli_abort(
      "{.field {label}_data} contains source{?s} not present elsewhere: {.val {unknown}}.",
      call = NULL
    )
  }

  if (!means_sds) {
    missing_cols <- setdiff(isotope_names, names(data))
    if (length(missing_cols) > 0) {
      cli::cli_abort(
        "Raw {.field {label}_data} is missing isotope column{?s}: {.field {missing_cols}}.",
        call = NULL
      )
    }
    Y <- as.matrix(data[, isotope_names, drop = FALSE])
    storage.mode(Y) <- "double"
    if (anyNA(Y)) {
      cli::cli_abort("Raw {.field {label}_data} isotope columns contain missing values.", call = NULL)
    }
    src_idx <- match(data[[source_col]], source_names)
    list(
      mode = "raw",
      n = nrow(Y),
      source_idx = src_idx,
      Y = Y
    )
  } else {
    mean_cols <- paste0(isotope_names, "_mean")
    sd_cols <- paste0(isotope_names, "_sd")
    missing_cols <- setdiff(c(mean_cols, sd_cols), names(data))
    if (length(missing_cols) > 0) {
      cli::cli_abort(c(
        "Summarised {.field {label}_data} is missing column{?s}: {.field {missing_cols}}.",
        "i" = "Expected one row per source, with columns {.field {mean_cols}} and {.field {sd_cols}}."
      ), call = NULL)
    }
    if (anyDuplicated(data[[source_col]]) > 0) {
      cli::cli_abort(
        "Summarised {.field {label}_data} must have exactly one row per source.",
        call = NULL
      )
    }
    if (!setequal(data[[source_col]], source_names)) {
      cli::cli_abort(
        "Summarised {.field {label}_data} must contain exactly the source{?s}: {.val {source_names}}.",
        call = NULL
      )
    }
    data <- data[match(source_names, data[[source_col]]), , drop = FALSE]
    mean_mat <- as.matrix(data[, mean_cols, drop = FALSE])
    sd_mat <- as.matrix(data[, sd_cols, drop = FALSE])
    storage.mode(mean_mat) <- storage.mode(sd_mat) <- "double"
    dimnames(mean_mat) <- dimnames(sd_mat) <- NULL
    if (anyNA(mean_mat) || anyNA(sd_mat)) {
      cli::cli_abort("Summarised {.field {label}_data} mean/sd columns contain missing values.", call = NULL)
    }
    if (any(sd_mat < 0)) {
      cli::cli_abort("Summarised {.field {label}_data} sd columns must be non-negative.", call = NULL)
    }
    list(
      mode = "summary",
      mean = mean_mat,
      sd = sd_mat
    )
  }
}

#' Derive the concentration-dependence matrix from `<isotope>_conc` columns
#' embedded in `source_data`, or `NULL` if `conc_dep` is `FALSE`.
#'
#' @param source_data Source data frame (raw or summarised, matching
#'   `source_means_sds`) -- the same one passed to `prep_iso_table()`. When
#'   `conc_dep = TRUE`, must additionally have a `<isotope>_conc` column for
#'   every isotope in `isotope_names`: one value per raw sample (averaged
#'   per source) if `source_means_sds = FALSE`, or one value per source
#'   (used directly) if `source_means_sds = TRUE`. Since each isotope's
#'   concentration is a proportion of a source's total mass, a source's
#'   values cannot sum to more than 1 across isotopes.
#' @param isotope_names Character vector of isotope names.
#' @param source_names Character vector of the full set of source names, in
#'   the order rows of the returned matrix should follow.
#' @param source_means_sds Logical; is `source_data` summarised (one row per
#'   source, `TRUE`) or raw replicate samples (`FALSE`)?
#' @param conc_dep Logical; enable concentration dependence (`TRUE`) by
#'   reading it from `source_data`'s `<isotope>_conc` columns, or disable it
#'   (`FALSE`).
#' @param source_col Name of the source-identifier column in `source_data`.
#' @return Numeric `length(source_names) x length(isotope_names)` matrix, or
#'   `NULL` if `conc_dep` is `FALSE`.
#' @noRd
prep_conc_dep <- function(source_data, isotope_names, source_names, source_means_sds,
                           conc_dep, source_col = "Source") {
  if (!isTRUE(conc_dep) && !isFALSE(conc_dep)) {
    cli::cli_abort(
      c(
        "{.arg conc_dep} must be {.code TRUE} or {.code FALSE}.",
        "i" = "Concentration values are read directly from {.field <isotope>_conc} column(s) in {.field source_data}, rather than a separate data frame."
      ),
      call = NULL
    )
  }
  if (!conc_dep) return(NULL)

  conc_cols <- paste0(isotope_names, "_conc")
  missing_cols <- setdiff(conc_cols, names(source_data))
  if (length(missing_cols) > 0) {
    cli::cli_abort(
      "{.code conc_dep = TRUE} requires concentration column{?s} in {.field source_data}: {.field {missing_cols}}.",
      call = NULL
    )
  }
  src_col_chr <- as.character(source_data[[source_col]])

  if (source_means_sds) {
    if (anyDuplicated(src_col_chr) > 0) {
      cli::cli_abort(
        "{.field source_data} must have exactly one row per source when {.code source_means_sds = TRUE}.",
        call = NULL
      )
    }
    m <- as.matrix(source_data[match(source_names, src_col_chr), conc_cols, drop = FALSE])
  } else {
    conc_mat <- as.matrix(source_data[, conc_cols, drop = FALSE])
    storage.mode(conc_mat) <- "double"
    # do.call(rbind, lapply(...)) rather than t(vapply(...)): vapply()
    # silently returns a plain vector, not a matrix, when FUN.VALUE has
    # length 1 (a single isotope), which t() would then turn into a 1 x K
    # row instead of the intended K x 1 column.
    m <- do.call(rbind, lapply(source_names, function(s) {
      colMeans(conc_mat[src_col_chr == s, , drop = FALSE])
    }))
  }

  storage.mode(m) <- "double"
  dimnames(m) <- NULL
  if (anyNA(m)) {
    cli::cli_abort("Concentration column(s) in {.field source_data} contain missing values.", call = NULL)
  }
  if (any(m <= 0 | m > 1)) {
    cli::cli_abort(
      "{.field source_data} concentration values must be in {.val (0, 1]} (elemental concentration proportions).",
      call = NULL
    )
  }
  row_sums <- rowSums(m)
  if (any(row_sums > 1)) {
    cli::cli_abort(
      c(
        "{.field source_data} concentration values sum to more than 1 across isotopes for source{?s}: {.val {source_names[row_sums > 1]}}.",
        "i" = "Each isotope's concentration is a proportion of that source's total mass, so they cannot sum to more than 1."
      ),
      call = NULL
    )
  }
  m
}
