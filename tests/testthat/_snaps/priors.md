# bsimms_prior errors on an invalid class

    Code
      bsimms_prior(prior = "normal(0, 1)", class = "banana")
    Condition
      Error in `bsimms_prior()`:
      ! `class` must be one of "b", "p_global", "sd", "cor", "sigma", "resid_prop", "source_mean", "source_sd", "tdf_mean", "tdf_sd", "source_cor", or "resid_cor", not "banana".

# bsimms_prior errors on a non-string prior

    Code
      bsimms_prior(prior = 1, class = "b")
    Condition
      Error:
      ! `prior` must be a single character string of Stan code, e.g. `"normal(0, 1)"`.

---

    Code
      bsimms_prior(prior = c("a", "b"), class = "b")
    Condition
      Error:
      ! `prior` must be a single character string of Stan code, e.g. `"normal(0, 1)"`.

# print.bsimms_prior prints a table with headers

    Code
      print(p)
    Output
                    prior class coef resp  group
             normal(0, 2)     b                 
       student_t(3, 0, 1)    sd           Region

# merge_bsimms_prior errors on a non-bsimms_prior user argument

    Code
      merge_bsimms_prior(default, data.frame(x = 1))
    Condition
      Error:
      ! `prior` must be built with `bsimms_prior()` (optionally combined with `c()`).

# merge_bsimms_prior errors on an unrecognised coef/resp/group

    Code
      merge_bsimms_prior(default, bsimms_prior("normal(0, 1)", class = "b", coef = "SexZZZ"),
      spec_re)
    Condition
      Error:
      ! `prior`: unrecognised coef "SexZZZ" for class "b"; must be one of "SexM".

---

    Code
      merge_bsimms_prior(default, bsimms_prior("2", class = "p_global", group = "Beavr"),
      spec_re)
    Condition
      Error:
      ! `prior`: unrecognised group "Beavr" for class "p_global"; must be one of "Beaver" and "Deer".

---

    Code
      merge_bsimms_prior(default, bsimms_prior("student_t(3, 0, 1)", class = "sd",
        group = "Regionn"), spec_re)
    Condition
      Error:
      ! `prior`: unrecognised group "Regionn" for class "sd"; must be one of "Region".

# merge_bsimms_prior errors clearly when a class has no valid terms in this model

    Code
      merge_bsimms_prior(default, bsimms_prior("student_t(3, 0, 1)", class = "sd",
        group = "Region"), spec)
    Condition
      Error:
      ! `prior`: class "sd" has no valid group in this model (no matching terms).

# select_prior errors when the class is entirely absent

    Code
      select_prior(df, "sd")
    Condition
      Error:
      ! No prior available for class "sd"; this should not happen — please report a bug.

# select_prior errors when no row matches any specificity level

    Code
      select_prior(df, "sd", coef = "y", group = "B")
    Condition
      Error:
      ! No matching prior found for class "sd", coef "y", resp "", group "B".

