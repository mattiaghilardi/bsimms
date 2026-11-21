# resolve_backend errors on an explicit backend that isn't installed

    Code
      resolve_backend("not_a_real_backend_pkg")
    Condition
      Error in `resolve_backend()`:
      ! The package "not_a_real_backend_pkg" is required to use `backend = "not_a_real_backend_pkg"`.

# resolve_backend errors when auto finds neither backend installed

    Code
      resolve_backend("auto")
    Condition
      Error:
      ! Fitting requires either cmdstanr or rstan to be installed.
      * `install.packages("rstan")`
      * or
      * `install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))`
      * `cmdstanr::install_cmdstan()`
      i You can still call `make_stancode()`/`make_standata()` and fit the model yourself.

# bsimm errors on an invalid error_structure or backend

    Code
      bsimm(formula = ~1, mixture_data = mixture_data, source_data = source_data,
        tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), error_structure = "banana",
        source_means_sds = TRUE)
    Condition
      Error in `bsimm()`:
      ! `error_structure` must be one of "process_residual", "process_only", or "residual_only", not "banana".

---

    Code
      bsimm(formula = ~1, mixture_data = mixture_data, source_data = source_data,
        tdf_data = tdf_data, isotope_names = c("d13C", "d15N"), backend = "banana",
        source_means_sds = TRUE)
    Condition
      Error in `bsimm()`:
      ! `backend` must be one of "auto", "cmdstanr", or "rstan", not "banana".

