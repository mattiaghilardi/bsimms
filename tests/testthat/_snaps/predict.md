# posterior_proportions errors when ndraws exceeds the number of draws available

    Code
      posterior_proportions(fit, ndraws = 1000)
    Condition
      Error:
      ! `ndraws` (1000) cannot exceed the number of posterior draws available (100).

# posterior_proportions errors by default (re_formula = NULL) when newdata is missing a group-level column

    Code
      posterior_proportions(fit, newdata = newdata)
    Condition
      Error:
      ! `newdata` is missing column Region, needed for its group-level term.

# posterior_proportions errors on unseen group levels in newdata

    Code
      posterior_proportions(fit, newdata = newdata)
    Condition
      Error:
      ! `newdata` contains levels of Region not seen when fitting the model: "C".
      i Set `allow_new_levels = TRUE` to predict for new levels, or use `re_formula` to exclude this term.

# posterior_proportions treats NA in a grouping column as a new level

    Code
      posterior_proportions(fit, newdata = newdata)
    Condition
      Error:
      ! `newdata` contains levels of Region not seen when fitting the model: NA.
      i Set `allow_new_levels = TRUE` to predict for new levels, or use `re_formula` to exclude this term.

# posterior_proportions errors when newdata is missing a fixed-effect column

    Code
      posterior_proportions(fit, newdata = newdata)
    Condition
      Error:
      ! `newdata` is missing column needed to build the design matrix: elevation.

# posterior_proportions errors when re_formula names a term not in the model

    Code
      posterior_proportions(fit, newdata = newdata, re_formula = ~ (1 | Foo))
    Condition
      Error:
      ! `re_formula` refers to group-level term not in the fitted model: Foo.

# posterior_proportions errors when re_formula is not NULL/NA/a formula

    Code
      posterior_proportions(fit, newdata = newdata, re_formula = "banana")
    Condition
      Error:
      ! `re_formula` must be `NULL`, `NA`, `~0`, or a formula naming group-level term(s), e.g. `~ (1 | Region)`.

# posterior_epred errors on an invalid resp

    Code
      rstantools::posterior_epred(fit, resp = "banana")
    Condition
      Error in `rstantools::posterior_epred()`:
      ! `resp` must be one of "d13C" or "d15N", not "banana".

# posterior_epred propagates posterior_proportions's newdata validation errors

    Code
      rstantools::posterior_epred(fit, newdata = newdata)
    Condition
      Error:
      ! `newdata` is missing column Region, needed for its group-level term.

# posterior_proportions errors by default when the inner nesting term's data is omitted

    Code
      posterior_proportions(nested_fit, newdata = data.frame(Site = factor("S1")))
    Condition
      Error:
      ! `newdata` is missing column needed for group-level term "Individual:Site":
      x Individual.

# posterior_proportions errors when the inner nesting level is supplied without the outer one

    Code
      posterior_proportions(nested_fit, newdata = data.frame(Individual = factor("I1")))
    Condition
      Error:
      ! `newdata` is missing column needed for group-level term "Individual:Site":
      x Site.

# posterior_proportions errors on an unseen Individual/Site combination, referencing the actual columns

    Code
      posterior_proportions(nested_fit, newdata = data.frame(Site = factor("S1"),
      Individual = factor("I3")))
    Condition
      Error:
      ! `newdata` contains a combination of Individual and Site not seen when fitting the model: "I3:S1".
      i Set `allow_new_levels = TRUE` to predict for new levels, or use `re_formula` to exclude this term.

# posterior_proportions errors on a non-factor newdata column for a factor fixed effect

    Code
      posterior_proportions(factor_fit, newdata = data.frame(Region = 1))
    Condition
      Error:
      ! Region must be a factor or character in `newdata`, as it was when fitting the model.

# posterior_proportions errors on an unseen factor level for a fixed effect

    Code
      posterior_proportions(factor_fit, newdata = data.frame(Region = factor("C",
        levels = c("A", "B", "C"))))
    Condition
      Error:
      ! `newdata` contains levels of Region not seen when fitting the model: "C".
      i Valid levels for Region are: "A" and "B".

