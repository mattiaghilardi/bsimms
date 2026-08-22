# bsimms_draws errors on a non-bsimms_fit object

    Code
      bsimms_draws(list())
    Condition
      Error:
      ! `object` must be a <bsimms_fit> object.

# extract_array_draws errors when the requested draws are missing

    Code
      extract_array_draws(dm, "p_global", dim1 = 5)
    Condition
      Error:
      ! Could not find draws for: p_global[3], p_global[4], and p_global[5].

