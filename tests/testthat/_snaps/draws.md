# bsimms_draws errors on a non-bsimms_fit object

    Code
      bsimms_draws(list())
    Condition
      Error:
      ! `object` must be a <bsimms_fit> object.

# draws_long errors when arr lacks dimnames on the 3rd dimension

    Code
      draws_long(arr)
    Condition
      Error:
      ! `arr` must be a `[n_draws, n_obs, n_var]` array with variable names attached as the 3rd dimension's dimnames.

# draws_long errors when arr is not 3-dimensional

    Code
      draws_long(matrix(1:4, 2, 2))
    Condition
      Error:
      ! `arr` must be a `[n_draws, n_obs, n_var]` array with variable names attached as the 3rd dimension's dimnames.

# draws_long errors when var_col and value_col are identical

    Code
      draws_long(arr, var_col = "x", value_col = "x")
    Condition
      Error:
      ! `var_col` and `value_col` must be different.

# extract_array_draws errors when the requested draws are missing

    Code
      extract_array_draws(dm, "p_global", dim1 = 5)
    Condition
      Error:
      ! Could not find draws for: p_global[3], p_global[4], and p_global[5].

