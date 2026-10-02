# Parse an `lme4`-style formula for the source-proportion model

`bsimms` formulas describe how the (ILR-transformed) source proportions
depend on mixture-level covariates, using the same syntax as `lme4`
(parsed here via the `reformulas` package, which now houses this
formula-processing machinery upstream of `lme4`/`glmmTMB`): fixed-effect
terms as usual (`~ Sex + Season`), and group-level ("random-effect")
terms in parentheses with a bar (`(1 | Region)`,
`(1 + Season | Individual)`, `(x || Region)` for an uncorrelated slope
and intercept). Formulas must not have a left-hand side: the response
(source proportions) is never observed directly, and is instead inferred
from the isotope mixture likelihood.

## Usage

``` r
parse_bsimms_formula(formula, data)
```

## Arguments

- formula:

  A one-sided formula, e.g. `~ Sex + Season + (1 | Region)`.

- data:

  The mixture data frame.

## Value

A list with elements `formula`, `fixed_formula`, `X` (fixed-effect
design matrix), `fixed_names`, `fixed_frame` (model frame of the
fixed-effect covariates, one column per named variable, character
columns coerced to factor), and `re_terms` (a list, one element per
group-level term, each with `group`, `term_names`, `Z`, `group_idx`,
`group_levels`).

## Details

All variables referenced in `formula` must be columns of the *mixture*
data set (source proportions describe mixture samples, and covariates
such as sex, age class, season, or capture site are properties of the
mixture sample, not the source).

There is no limit on the number or type of fixed-effect covariates
(continuous, factor, interactions, ...) and no limit on the number of
independent, crossed, or nested group-level (hierarchical) terms:
anything `lme4::lmer()` could parse on the right-hand side of a formula
is supported here, including `/` for nested terms (e.g.
`(1 | site/individual)`, expanded exactly as in `lme4` into
`(1 | site) + (1 | individual:site)`) and `:` for crossed terms (e.g.
`(1 | site:individual)`).

## Examples

``` r
d <- data.frame(Sex = c("M", "F", "F", "M"), Region = c("A", "A", "B", "B"))
pf <- parse_bsimms_formula(~ Sex + (1 | Region), data = d)
pf$fixed_names
#> [1] "(Intercept)" "SexM"       
pf$re_terms[[1]]$group
#> [1] "Region"
```
