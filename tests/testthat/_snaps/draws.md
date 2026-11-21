# bsimms_draws errors on a non-bsimms_fit object

    Code
      bsimms_draws(list())
    Condition
      Error:
      ! `object` must be a <bsimms_fit> object.

# extract_array_draws reshapes a matrix parameter

    Code
      extract_array_draws(dm, "source_mean", dim1 = 2, dim2 = 5)
    Condition
      Error:
      ! Could not find draws for: source_mean[1,3], source_mean[1,4], and source_mean[1,5].

# extract_array_of_matrices reshapes an array of matrices and errors when draws are missing

    Code
      extract_array_of_matrices(dm, "Lcorr_source", dim1 = 5, dim2 = 2)
    Condition
      Error:
      ! Could not find draws for: Lcorr_source[3,1,1] and Lcorr_source[3,1,2].

# subset_ndraws errors on an invalid ndraws

    Code
      subset_ndraws(dm, ndraws = -1)
    Condition
      Error:
      ! `ndraws` must be a single positive integer.

---

    Code
      subset_ndraws(dm, ndraws = "a")
    Condition
      Error:
      ! `ndraws` must be a single positive integer.

---

    Code
      subset_ndraws(dm, ndraws = 1.5)
    Condition
      Error:
      ! `ndraws` must be a single positive integer.

---

    Code
      subset_ndraws(dm, ndraws = c(1, 2))
    Condition
      Error:
      ! `ndraws` must be a single positive integer.

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

