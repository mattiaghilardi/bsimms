# bsimm errors on an invalid error_structure or backend

    Code
      bsimm(formula = ~1, mixture_data = mixture_data, source_data = source_data,
        tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), error_structure = "banana")
    Condition
      Error in `bsimm()`:
      ! `error_structure` must be one of "process_residual", "process_only", or "residual_only", not "banana".

---

    Code
      bsimm(formula = ~1, mixture_data = mixture_data, source_data = source_data,
        tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), backend = "banana")
    Condition
      Error in `bsimm()`:
      ! `backend` must be one of "auto", "cmdstanr", or "rstan", not "banana".

