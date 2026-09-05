#' Parse an `lme4`-style formula for the source-proportion model
#'
#' `bsimms` formulas describe how the (ILR-transformed) source proportions
#' depend on mixture-level covariates, using the same syntax as `lme4`
#' (parsed here via the `reformulas` package, which now houses this
#' formula-processing machinery upstream of `lme4`/`glmmTMB`): fixed-effect
#' terms as usual (`~ Sex + Season`), and group-level ("random-effect")
#' terms in parentheses with a bar (`(1 | Region)`, `(1 + Season |
#' Individual)`, `(x || Region)` for an uncorrelated slope and intercept).
#' Formulas must not have a left-hand side: the response (source
#' proportions) is never observed directly, and is instead inferred from
#' the isotope mixture likelihood.
#'
#' All variables referenced in `formula` must be columns of the *mixture*
#' data set (source proportions describe mixture samples, and covariates
#' such as sex, age class, season, or capture site are properties of the
#' mixture sample, not the source).
#'
#' There is no limit on the number or type of fixed-effect covariates
#' (continuous, factor, interactions, ...) and no limit on
#' the number of independent, crossed, or nested group-level (hierarchical)
#' terms: anything `lme4::lmer()` could parse on the right-hand side of a
#' formula is supported here, including `/` for nested terms (e.g.
#' `(1 | site/individual)`, expanded exactly as in `lme4` into
#' `(1 | site) + (1 | individual:site)`) and `:` for crossed terms (e.g.
#' `(1 | site:individual)`).
#'
#' @param formula A one-sided formula, e.g. `~ Sex + Season + (1 | Region)`.
#' @param data The mixture data frame.
#' @return A list with elements `formula`, `fixed_formula`, `X` (fixed-effect
#'   design matrix), `fixed_names`, `fixed_frame` (model frame of the
#'   fixed-effect covariates, one column per named variable, character
#'   columns coerced to factor), and `re_terms` (a list, one element per
#'   group-level term, each with `group`, `term_names`, `Z`, `group_idx`,
#'   `group_levels`).
#' @export
#' @examples
#' d <- data.frame(Sex = c("M", "F", "F", "M"), Region = c("A", "A", "B", "B"))
#' pf <- parse_bsimms_formula(~ Sex + (1 | Region), data = d)
#' pf$fixed_names
#' pf$re_terms[[1]]$group
parse_bsimms_formula <- function(formula, data) {
  if (!inherits(formula, "formula")) {
    cli::cli_abort(
      "{.arg formula} must be a formula, e.g. {.code ~ Sex + (1 | Region)}.",
      call = NULL
    )
  }
  if (length(formula) == 3) {
    cli::cli_abort(
      c(
        paste0(
          "{.arg formula} must not have a left-hand side (e.g. use ",
          "{.code ~ Sex}, not {.code p ~ Sex})."
        ),
        "i" = paste0(
          "Source proportions are inferred from the isotope mixture, ",
          "not supplied directly."
        )
      ),
      call = NULL
    )
  }
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame.", call = NULL)
  }

  # Expand `(x || g)` into independent single-bar terms, exactly as lme4/
  # glmmTMB do internally; after this every remaining bar term is a plain
  # `term | group`.
  formula <- reformulas::expandDoubleVerts(formula)

  bars <- reformulas::findbars(formula)
  fixed_formula <- reformulas::nobars(formula)
  # nobars() on a one-sided formula can return `~1`; keep an explicit
  # intercept-only design in that case.
  X <- stats::model.matrix(fixed_formula, data = data)
  fixed_names <- colnames(X)
  fixed_frame <- stats::model.frame(fixed_formula, data = data)
  fixed_frame[] <- lapply(fixed_frame, function(col) {
    if (is.character(col)) factor(col) else col
  })

  # Character columns must be coerced to factor before evaluating a
  # group-level term's grouping expression: base R's `:` only computes the
  # factor interaction (e.g. `a:b` -> "x:p", "x:q", ...) when both operands
  # are already factors -- on character vectors it silently falls back to
  # arithmetic sequence generation instead, so a crossed grouping factor
  # like `(1 | a:b)` would otherwise error or misbehave whenever `a`/`b`
  # are stored as character (the common case, since R defaults to
  # `stringsAsFactors = FALSE`).
  data_for_groups <- data
  data_for_groups[] <- lapply(data_for_groups, function(col) {
    if (is.character(col)) factor(col) else col
  })

  re_terms <- list()
  if (!is.null(bars)) {
    for (b in bars) {
      b_str <- rlang::expr_deparse(b)
      parts <- strsplit(b_str, "\\|")[[1]]
      if (length(parts) != 2) {
        cli::cli_abort(
          "Could not parse group-level term {.code {b_str}}.",
          call = NULL
        )
      }
      term_str <- trimws(parts[1])
      group_str <- trimws(parts[2])
      term_formula <- stats::as.formula(paste("~", term_str))

      group_val <- rlang::eval_tidy(
        rlang::parse_expr(group_str),
        data = data_for_groups,
        env = parent.frame()
      )
      group_factor <- factor(group_val)
      if (any(is.na(group_factor))) {
        cli::cli_abort(
          "Grouping factor {.field {group_str}} contains missing values.",
          call = NULL
        )
      }

      Zg <- stats::model.matrix(term_formula, data = data)

      re_terms[[length(re_terms) + 1]] <- list(
        group = group_str,
        term_formula = term_formula,
        term_names = colnames(Zg),
        Z = Zg,
        group_idx = as.integer(group_factor),
        group_levels = levels(group_factor)
      )
    }
  }

  # Give each group-level term a unique, valid-Stan-identifier label, even
  # when the same grouping factor is used in more than one term (e.g. from
  # an expanded `||`).
  if (length(re_terms) > 0) {
    group_names <- vapply(re_terms, function(x) x$group, character(1))
    dup <- stats::ave(seq_along(group_names), group_names, FUN = seq_along)
    labels <- ifelse(
      dup == 1,
      make_stan_name(group_names),
      paste0(make_stan_name(group_names), "_", dup)
    )
    for (i in seq_along(re_terms)) {
      re_terms[[i]]$label <- labels[i]
    }
  }

  list(
    formula = formula,
    fixed_formula = fixed_formula,
    X = X,
    fixed_names = fixed_names,
    fixed_frame = fixed_frame,
    re_terms = re_terms
  )
}

#' Sanitise a character vector into valid Stan identifiers: non-alphanumeric
#' characters (besides `_`) become `_`, and any result not starting with a
#' letter is prefixed with `g_`. Used to turn a group-level term's grouping
#' factor name (e.g. `"Site:Individual"`) into a name suitable for a Stan
#' parameter/data variable (e.g. `sd_re_<label>`).
#'
#' @param x Character vector.
#' @return Character vector of the same length, valid as Stan identifiers.
#' @noRd
make_stan_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- ifelse(grepl("^[A-Za-z]", x), x, paste0("g_", x))
  x
}
