# infer_source_names errors on a non-string source_col

    Code
      infer_source_names(source_data = sources, tdf_data = NULL, source_col = c(
        "Source", "Other"))
    Condition
      Error:
      ! `source_col` must be a single string.

---

    Code
      infer_source_names(source_data = sources, tdf_data = NULL, source_col = NULL)
    Condition
      Error:
      ! `source_col` must be a single string.

# infer_source_names errors on missing source_col

    Code
      infer_source_names(source_data = data.frame(x = 1), tdf_data = NULL)
    Condition
      Error:
      ! Source data must have a Source column.

---

    Code
      infer_source_names(source_data = data.frame(Source = "A"), tdf_data = data.frame(
        x = 1))
    Condition
      Error:
      ! TDF data must have a Source column.

# infer_source_names errors on source/tdf name mismatch

    Code
      infer_source_names(source_data = sources, tdf_data = tdf)
    Condition
      Error:
      ! Source names in source_data and tdf_data do not match.
      x source_data: "Beaver" and "Deer"
      x tdf_data: "Beaver" and "Hare"

# prep_mixture_isotopes errors on missing isotope columns

    Code
      prep_mixture_isotopes(mixture_data = data.frame(d13C = 1), isotope_names = c(
        "d13C", "d15N"))
    Condition
      Error:
      ! Mixture data is missing isotope column: d15N.

# prep_mixture_isotopes errors on missing values

    Code
      prep_mixture_isotopes(mixture_data = d, isotope_names = c("d13C", "d15N"))
    Condition
      Error:
      ! Mixture isotope data contain missing values; remove or impute them first.

# prep_iso_table errors on an invalid label

    Code
      prep_iso_table(data = d, isotope_names = "d13C", source_names = "Beaver",
        means_sds = FALSE, label = "banana")
    Condition
      Error in `prep_iso_table()`:
      ! `label` must be one of "source" or "tdf", not "banana".

# prep_iso_table (raw) errors on missing source_col

    Code
      prep_iso_table(data = data.frame(x = 1), isotope_names = "d13C", source_names = "Beaver",
      means_sds = FALSE)
    Condition
      Error:
      ! source_data must have a Source column.

# prep_iso_table (raw) errors on unknown source

    Code
      prep_iso_table(data = d, isotope_names = "d13C", source_names = c("Beaver",
        "Deer"), means_sds = FALSE)
    Condition
      Error:
      ! source_data contains source not present elsewhere: "Fox".

# prep_iso_table (raw) errors on missing isotope columns

    Code
      prep_iso_table(data = d, isotope_names = c("d13C", "d15N"), source_names = "Beaver",
      means_sds = FALSE)
    Condition
      Error:
      ! Raw source_data is missing isotope column: d15N.

# prep_iso_table (raw) errors on missing values

    Code
      prep_iso_table(data = d, isotope_names = "d13C", source_names = "Beaver",
        means_sds = FALSE)
    Condition
      Error:
      ! Raw source_data isotope columns contain missing values.

# prep_iso_table (summary) errors on missing mean/sd columns

    Code
      prep_iso_table(data = d, isotope_names = "d13C", source_names = "Beaver",
        means_sds = TRUE)
    Condition
      Error:
      ! Summarised source_data is missing column: d13C_sd.
      i Expected one row per source, with columns d13C_mean and d13C_sd.

# prep_iso_table (summary) errors on duplicated source rows

    Code
      prep_iso_table(data = d, isotope_names = "d13C", source_names = "Beaver",
        means_sds = TRUE)
    Condition
      Error:
      ! Summarised source_data must have exactly one row per source.

# prep_iso_table (summary) errors when source set doesn't match exactly

    Code
      prep_iso_table(data = d, isotope_names = "d13C", source_names = c("Beaver",
        "Deer"), means_sds = TRUE)
    Condition
      Error:
      ! Summarised source_data must contain exactly the source: "Beaver" and "Deer".

# prep_iso_table (summary) errors on missing values

    Code
      prep_iso_table(data = d, isotope_names = "d13C", source_names = "Beaver",
        means_sds = TRUE)
    Condition
      Error:
      ! Summarised source_data mean/sd columns contain missing values.

# prep_iso_table (summary) errors on negative sd

    Code
      prep_iso_table(data = d, isotope_names = "d13C", source_names = "Beaver",
        means_sds = TRUE)
    Condition
      Error:
      ! Summarised source_data sd columns must be non-negative.

# prep_conc_dep errors when conc_dep is not TRUE/FALSE

    Code
      prep_conc_dep(source_data = data.frame(Source = "Beaver"), isotope_names = "d13C",
      source_names = "Beaver", source_means_sds = FALSE, conc_dep = "yes")
    Condition
      Error:
      ! `conc_dep` must be `TRUE` or `FALSE`.
      i Concentration values are read directly from <isotope>_conc column(s) in source_data, rather than a separate data frame.

# prep_conc_dep errors on missing concentration columns

    Code
      prep_conc_dep(source_data = d, isotope_names = "d13C", source_names = "Beaver",
        source_means_sds = FALSE, conc_dep = TRUE)
    Condition
      Error:
      ! `conc_dep = TRUE` requires concentration column in source_data: d13C_conc.

# prep_conc_dep (summary) errors on duplicated source rows

    Code
      prep_conc_dep(source_data = d, isotope_names = "d13C", source_names = "Beaver",
        source_means_sds = TRUE, conc_dep = TRUE)
    Condition
      Error:
      ! source_data must have exactly one row per source when `source_means_sds = TRUE`.

# prep_conc_dep errors on missing concentration values

    Code
      prep_conc_dep(source_data = d, isotope_names = "d13C", source_names = "Beaver",
        source_means_sds = FALSE, conc_dep = TRUE)
    Condition
      Error:
      ! Concentration column(s) in source_data contain missing values.

# prep_conc_dep errors on values outside (0, 1]

    Code
      prep_conc_dep(source_data = d, isotope_names = "d13C", source_names = "Beaver",
        source_means_sds = TRUE, conc_dep = TRUE)
    Condition
      Error:
      ! source_data concentration values must be in "(0, 1]" (elemental concentration proportions).

# prep_conc_dep errors when a source's concentrations sum to more than 1

    Code
      prep_conc_dep(source_data = d, isotope_names = c("d13C", "d15N"), source_names = "Beaver",
      source_means_sds = TRUE, conc_dep = TRUE)
    Condition
      Error:
      ! source_data concentration values sum to more than 1 across isotopes for source: "Beaver".
      i Each isotope's concentration is a proportion of that source's total mass, so they cannot sum to more than 1.

