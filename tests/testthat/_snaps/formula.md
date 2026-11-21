# a two-sided formula is rejected

    Code
      parse_bsimms_formula(formula = y ~ x, data = d)
    Condition
      Error:
      ! `formula` must not have a left-hand side (e.g. use `~ Sex`, not `p ~ Sex`).
      i Source proportions are inferred from the isotope mixture, not supplied directly.

# a non-formula is rejected

    Code
      parse_bsimms_formula(formula = "~x", data = d)
    Condition
      Error:
      ! `formula` must be a formula, e.g. `~ Sex + (1 | Region)`.

# non-data-frame data is rejected

    Code
      parse_bsimms_formula(formula = ~x, data = list(x = 1:3))
    Condition
      Error:
      ! `data` must be a data frame.

# a chained bar term is rejected

    Code
      parse_bsimms_formula(formula = ~ 1 + (1 | a | b), data = d)
    Condition
      Error:
      ! Could not parse group-level term `1 | a | b`.

# a grouping factor with missing values is rejected

    Code
      parse_bsimms_formula(formula = ~ 1 + (1 | Region), data = d)
    Condition
      Error:
      ! Grouping factor Region contains missing values.

