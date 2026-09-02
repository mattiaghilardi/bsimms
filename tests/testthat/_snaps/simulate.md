# unknown n_levels/n_groups names are rejected

    Code
      simulate_bsimms_data(~Sex, n_mixture_obs = 10, n_levels = list(banana = 2))
    Condition
      Error:
      ! `n_levels` names variable not in `formula`: "banana".

---

    Code
      simulate_bsimms_data(~ (1 | Region), n_mixture_obs = 10, n_groups = list(
        banana = 2))
    Condition
      Error:
      ! `n_groups` names variable not in `formula`: "banana".

# a grouping factor missing from n_groups is rejected

    Code
      simulate_bsimms_data(~ (1 | Region), n_mixture_obs = 10)
    Condition
      Error:
      ! `n_groups` must include an entry for every grouping factor in `formula`: missing "Region".

# random slopes are rejected

    Code
      simulate_bsimms_data(~ (elevation | Region), n_mixture_obs = 10, n_groups = list(
        Region = 2))
    Condition
      Error:
      ! `formula` contains a random slope (`elevation | Region`).
      i `simulate_bsimms_data()` only supports random intercepts, e.g. `(1 | Group)`.

# n_levels/n_groups entries must be integers >= 2

    Code
      simulate_bsimms_data(~Sex, n_mixture_obs = 10, n_levels = list(Sex = 1))
    Condition
      Error:
      ! `n_levels` entries must be single integers >= 2.

# too many levels for the available observations is rejected

    Code
      simulate_bsimms_data(~Region, n_mixture_obs = 2, n_levels = list(Region = 5))
    Condition
      Error:
      ! Cannot split 2 observations across 5 levels: each level needs at least one observation.

# p_global can be overridden and is validated

    Code
      simulate_bsimms_data(n_mixture_obs = 10, n_sources = 3, p_global = c(0.2, 0.5,
        0.5))
    Condition
      Error:
      ! `p_global` must be a vector of 3 strictly positive proportions summing to 1.

# n_source_obs/n_tdf_obs reject the wrong length or non-positive values

    Code
      simulate_bsimms_data(n_mixture_obs = 10, n_sources = 3, n_source_obs = c(5, 20))
    Condition
      Error:
      ! `n_source_obs` must be a single positive integer or a positive integer vector of length 3 (n_sources).

---

    Code
      simulate_bsimms_data(n_mixture_obs = 10, n_sources = 3, n_source_obs = c(5, -1,
        8))
    Condition
      Error:
      ! `n_source_obs` must be a single positive integer or a positive integer vector of length 3 (n_sources).

# an underdetermined system (n_sources > n_isotopes + 1) warns

    Code
      invisible(simulate_bsimms_data(n_mixture_obs = 10, n_sources = 4, n_isotopes = 2,
        seed = 1))
    Condition
      Warning:
      4 sources with only 2 isotopes is an underdetermined mixing system (`n_sources > n_isotopes + 1`).
      i Source proportions will rely more heavily on the prior than the data.

