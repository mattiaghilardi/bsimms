# plot.bsimms_fit errors when no parameters match

    Code
      plot(fit, variable = "banana", plot = FALSE)
    Condition
      Error:
      ! Can't find the following variable(s) in the output: banana

# pp_check.bsimms_fit requires resp for multi-isotope models

    Code
      bayesplot::pp_check(fit)
    Condition
      Error:
      ! Model has multiple isotopes ("d13C" and "d15N"); specify `resp` to select one.

# pp_check.bsimms_fit errors on an invalid type

    Code
      bayesplot::pp_check(fit, resp = "d13C", type = "banana")
    Condition
      Error:
      ! `type` must be one of "bars", "bars_grouped", "boxplot", "calibration", "calibration_grouped", "calibration_overlay", "calibration_overlay_grouped", "dens", "dens_overlay", "dens_overlay_grouped", "dots", "ecdf_overlay", "ecdf_overlay_grouped", "error_binned", "error_hist", "error_hist_grouped", "error_scatter", "error_scatter_avg", ..., "stat_grouped", and "violin_grouped".

# pp_check.bsimms_fit errors when ndraws exceeds the number of draws available

    Code
      bayesplot::pp_check(fit_1iso, ndraws = 1000)
    Condition
      Error:
      ! `ndraws` (1000) cannot exceed the number of posterior draws available (100).

