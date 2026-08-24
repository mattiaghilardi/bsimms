# bayes_R2 errors on an invalid resp

    Code
      rstantools::bayes_R2(fit, resp = "banana")
    Condition
      Error in `rstantools::bayes_R2()`:
      ! `resp` must be one of "d13C" or "d15N", not "banana".

