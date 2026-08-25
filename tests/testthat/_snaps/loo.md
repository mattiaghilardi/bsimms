# loo_compare errors when a precomputed object doesn't match criterion

    Code
      loo::loo_compare(fit2, l1, criterion = "waic")
    Condition
      Error:
      ! Cannot compare models evaluated with different criteria: "l1" must be <waic> object to match `criterion = 'waic'`.

