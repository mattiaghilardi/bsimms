# conditional_effects errors on a non-bsimms_fit object

    Code
      conditional_effects(list())
    Condition
      Error:
      ! `object` must be a <bsimms_fit> object.

# conditional_effects errors when the model has no fixed-effect covariates

    Code
      conditional_effects(fit_no_fixed)
    Condition
      Error:
      ! Model has no fixed-effect covariates to condition on.

# conditional_effects errors when effects names an unknown covariate

    Code
      conditional_effects(fit, effects = "banana")
    Condition
      Error:
      ! `effects` entry "banana" refers to covariate not in the model: "banana".

# ref_conditions errors on non-numeric values for a numeric covariate

    Code
      conditional_effects(fit, effects = "Region", ref_conditions = list(elevation = "high"))
    Condition
      Error:
      ! `ref_conditions` for elevation must be numeric, since it is a numeric covariate.

# ref_conditions errors on an unseen level for a factor covariate

    Code
      conditional_effects(fit, effects = "elevation", ref_conditions = list(Region = "X"))
    Condition
      Error:
      ! `ref_conditions` for Region contains a level not in the model: "X".
      i Valid levels for Region are: "A" and "B".

# ref_conditions errors on more than one value for a covariate

    Code
      conditional_effects(fit, effects = "Region", ref_conditions = list(elevation = c(
        100, 200)))
    Condition
      Error:
      ! `ref_conditions` for elevation must be a single value.

# ref_conditions errors on a misspelled covariate name instead of falling back to the default

    Code
      conditional_effects(fit, effects = "Region", ref_conditions = list(elevaton = 500))
    Condition
      Error:
      ! `ref_conditions` names covariate not in the model: "elevaton".

# int_conditions errors on a misspelled covariate name instead of falling back to the default

    Code
      conditional_effects(fit, effects = "elevation:Region", int_conditions = list(
        Regoin = "A"))
    Condition
      Error:
      ! `int_conditions` names covariate not in the model: "Regoin".

# conditional_effects's ... is forwarded to posterior_proportions

    Code
      conditional_effects(fit, effects = "elevation", ndraws = 1000)
    Condition
      Error:
      ! `ndraws` (1000) cannot exceed the number of posterior draws available (100).

# conditional_effects errors when re_formula names a term absent from the grid

    Code
      conditional_effects(fit, effects = "elevation", re_formula = ~ (1 | Pack))
    Condition
      Error:
      ! `newdata` is missing column Pack, needed for its group-level term.

# conditional_effects errors when an interaction names more than two covariates

    Code
      conditional_effects(fit, effects = "elevation:Region:Sex")
    Condition
      Error:
      ! `effects` entry "elevation:Region:Sex" names more than two covariates; only two-way interactions are supported.

# conditional_effects errors when an interaction names the same covariate twice

    Code
      conditional_effects(fit, effects = "elevation:elevation")
    Condition
      Error:
      ! `effects` entry "elevation:elevation" names the same covariate twice.

# conditional_effects errors when an interaction names an unknown covariate

    Code
      conditional_effects(fit, effects = "elevation:banana")
    Condition
      Error:
      ! `effects` entry "elevation:banana" refers to covariate not in the model: "banana".

# int_conditions errors on an unseen level for a factor moderator

    Code
      conditional_effects(fit, effects = "Region:Sex", int_conditions = list(Sex = "X"))
    Condition
      Error:
      ! `int_conditions` for Sex contains level not in the model: "X".
      i Valid levels for Sex are: "F" and "M".

# int_conditions errors on non-numeric values for a numeric moderator

    Code
      conditional_effects(fit, effects = "elevation:rainfall", int_conditions = list(
        rainfall = "high"))
    Condition
      Error:
      ! `int_conditions` for rainfall must be numeric, since it is a numeric covariate.

# conditional_effects errors on an invalid method

    Code
      conditional_effects(fit, effects = "elevation", method = "banana")
    Condition
      Error in `conditional_effects()`:
      ! `method` must be one of "posterior_proportions", "posterior_epred", or "posterior_predict", not "banana".

# conditional_effects requires resp for a multi-isotope model with method = posterior_epred

    Code
      conditional_effects(fit, effects = "elevation", method = "posterior_epred")
    Condition
      Error:
      ! Model has multiple isotopes ("d13C" and "d15N"); specify `resp` to select one.

# conditional_effects errors on an unknown resp

    Code
      conditional_effects(fit, effects = "elevation", method = "posterior_epred",
        resp = "banana")
    Condition
      Error in `conditional_effects()`:
      ! `resp` must be one of "d13C" or "d15N", not "banana".

