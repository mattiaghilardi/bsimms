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
      ! `type` must be one of "bars", "bars_grouped", "boxplot", "dens", "dens_overlay", "dens_overlay_grouped", "dots", "ecdf_overlay", "ecdf_overlay_grouped", "error_binned", "error_hist", "error_hist_grouped", "error_scatter", "error_scatter_avg", "error_scatter_avg_grouped", "error_scatter_avg_vs_x", "freqpoly", "freqpoly_grouped", ..., "stat_grouped", and "violin_grouped".

