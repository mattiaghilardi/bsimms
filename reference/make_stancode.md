# Generate the Stan code for a `bsimms` model

Builds a complete, human-readable Stan program implementing the stable
isotope mixing model described by `formula` and the supplied data:
nothing is hidden inside compiled internals, the generated `.stan` text
is the model, and it can be inspected, hand-edited, and compiled with
`cmdstanr` or `rstan` directly.

## Usage

``` r
make_stancode(
  formula,
  mixture_data,
  source_data,
  tdf_data,
  isotope_names,
  source_means_sds = FALSE,
  tdf_means_sds = TRUE,
  conc_dep = FALSE,
  error_structure = c("process_residual", "process_only", "residual_only"),
  prior = NULL,
  source_col = "Source"
)
```

## Arguments

- formula:

  A one-sided `lme4`-style formula for the source proportions, e.g.
  `~ Sex + Season + (1 + Season | Region)`. All variables referenced
  must be columns of `mixture_data`. The population-level baseline is
  not part of `formula`: it is `p_global`, a `simplex[K]` parameter with
  a `Dirichlet(alpha)` prior (as in `MixSIAR`), and `formula`'s terms
  are covariate deviations from it, added in ILR space (see
  [`bsimms_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_prior.md)'s
  `"p_global"` class to set `alpha`). A formula with no covariates (e.g.
  `~ 1`) is a valid Dirichlet-baseline-only model.

- mixture_data:

  Data frame of mixture observations (e.g. consumers, in a diet-mixing
  study): one row per individual (or per sample), with a column for each
  entry of `isotope_names`, plus any covariate/grouping columns used in
  `formula`.

- source_data:

  Source isotope data. If `source_means_sds = FALSE` (default), long
  format with one row per raw sample and columns `Source`,
  `isotope_names`. If `source_means_sds = TRUE`, one row per source with
  columns `Source`, `<isotope>_mean`, `<isotope>_sd` for every isotope.

- tdf_data:

  Trophic discrimination factor (diet-tissue discrimination) data, same
  layout convention as `source_data`, controlled by `tdf_means_sds`
  (default `TRUE`, since TDFs are usually taken from the literature as
  means/SDs).

- isotope_names:

  Character vector of isotope column names shared by `mixture_data`,
  `source_data` and `tdf_data`, e.g. `c("d13C", "d15N")`.

- source_means_sds:

  Logical; is `source_data` supplied as means/SDs (`TRUE`) or raw
  replicate samples (`FALSE`, default)? When raw, source means and SDs
  are estimated as part of the model, and their uncertainty is fully
  propagated into the posterior.

- tdf_means_sds:

  Logical; is `tdf_data` supplied as means/SDs (`TRUE`, default) or raw
  replicate samples (`FALSE`)? When raw, TDF means and SDs are estimated
  as part of the model and their uncertainty is fully propagated into
  the posterior. Unlike raw source data, though, no cross-isotope
  correlation is estimated for raw TDF data (there is no `"tdf_cor"`
  class in
  [`bsimms_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_prior.md)):
  the mixture likelihood only takes the multivariate
  (`multi_normal_cholesky`) form when *source* data is raw and there are
  2+ isotopes. This asymmetry is intentional (raw TDF data with enough
  replicates per source to usefully estimate a correlation is a rare
  combination in practice), not an oversight, but could be added
  symmetrically to `source_cor` if needed.

- conc_dep:

  Logical; enable elemental concentration dependence (default `FALSE`)?
  If `TRUE`, `source_data` must have a `<isotope>_conc` column for every
  isotope in `isotope_names`, giving that source's proportional
  elemental concentration (in `(0, 1]`) for that isotope's element – one
  value per raw sample (averaged per source) if
  `source_means_sds = FALSE`, or one value per source if
  `source_means_sds = TRUE`. Since these are proportions of the same
  source's total mass, they cannot sum to more than 1 across isotopes.

- error_structure:

  One of `"process_residual"` (default: source/TDF variance propagated
  into the mixture and scaled by an estimated multiplicative
  residual-error factor, i.e. the "Residual \* Process" error structure
  of Stock & Semmens 2016), `"process_only"` (only source/TDF variance
  propagated, no separate residual term, as in the original MixSIR
  model, Moore & Semmens 2008), or `"residual_only"` (source/TDF
  variance is not propagated into the mixture at all; all unexplained
  variance instead goes into an isotope-specific residual error term, as
  in the original SIAR model, Parnell et al. 2010).

- prior:

  Optional `bsimms_prior` object (see
  [`bsimms_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_prior.md))
  with one or more rows overriding the default priors. Unspecified
  parameters keep their (weakly informative, partly data-scaled)
  defaults; see
  [`bsimms_get_prior()`](https://mattiaghilardi.github.io/bsimms/reference/bsimms_get_prior.md).

- source_col:

  Name of the source-identifier column shared by `source_data` and
  `tdf_data`. Default `"Source"`.

## Value

A single character string of Stan code (class `bsimms_stancode`;
[`print()`](https://rdrr.io/r/base/print.html)s as plain text).

## References

Moore, J.W., & Semmens, B.X. (2008). Incorporating uncertainty and prior
information into stable isotope mixing models. *Ecology Letters*, 11(5),
470-480.
[doi:10.1111/j.1461-0248.2008.01163.x](https://doi.org/10.1111/j.1461-0248.2008.01163.x)

Parnell, A.C., Inger, R., Bearhop, S., & Jackson, A.L. (2010). Source
partitioning using stable isotopes: coping with too much variation.
*PLoS ONE*, 5(3), e9672.
[doi:10.1371/journal.pone.0009672](https://doi.org/10.1371/journal.pone.0009672)

Stock, B.C., & Semmens, B.X. (2016). Unifying error structures in
commonly used biotracer mixing models. *Ecology*, 97(10), 2562-2569.
[doi:10.1002/ecy.1517](https://doi.org/10.1002/ecy.1517)

## Examples

``` r
sim <- simulate_bsimms_data(
  ~ 1 + (1 | Region),
  n_mixture_obs = 20,
  source_names = c("Beaver", "Deer", "Hare"),
  isotope_names = c("d13C", "d15N"),
  n_groups = list(Region = 2),
  seed = 1
)
code <- make_stancode(
  sim$formula,
  mixture_data = sim$mixture_data,
  source_data = sim$source_data,
  tdf_data = sim$tdf_data,
  isotope_names = sim$isotope_names,
  source_means_sds = sim$source_means_sds,
  tdf_means_sds = sim$tdf_means_sds,
  conc_dep = sim$conc_dep,
  error_structure = sim$error_structure,
  source_col = sim$source_col
)
cat(code)
#> // Stan program generated by bsimms 0.0.0.9000
#> // sources: Beaver, Deer, Hare
#> // isotopes: d13C, d15N
#> // global proportions: p_global ~ Dirichlet(alpha), covariate effects added in ILR space
#> // fixed effects (covariate slopes, no intercept): (none)
#> // group-level terms: ((Intercept) | Region)
#> // error structure: process_residual
#> // source data: raw
#> // tdf data: summary
#> 
#> functions {
#>   /* Geometric mean of the first k parts of a composition
#>    * Used for the forward ILR transform of p_global (Egozcue et al. 2003,
#>    * eq. 25)
#>    * Args:
#>    *   p: a composition (simplex-scale vector)
#>    *   k: number of leading parts of p to average
#>    * Returns:
#>    *   the geometric mean of p[1:k]
#>    */
#>   real gmean(vector p, int k) {
#>     return exp(mean(log(p[1:k])));
#>   }
#> }
#> data {
#>   int<lower=1> N;  // number of mixture observations
#>   int<lower=1> J;  // number of isotopes/tracers
#>   int<lower=2> K;  // number of sources
#>   int<lower=1> D;  // K - 1, ILR dimension
#>   matrix[K, D] V;  // ILR basis matrix
#>   matrix[N, J] y;  // mixture isotope data
#>   vector<lower=0>[K] alpha_dirichlet;  // Dirichlet concentration for p_global
#>   int<lower=1> N_re_Region;  // levels of grouping factor 'Region'
#>   int<lower=1> M_re_Region;  // number of group-level terms for 'Region'
#>   matrix[N, M_re_Region] Z_re_Region;  // group-level design matrix for 'Region'
#>   array[N] int<lower=1, upper=N_re_Region> grp_re_Region;  // level index per observation
#>   int<lower=1> N_source_raw;  // number of raw source samples
#>   array[N_source_raw] int<lower=1, upper=K> source_idx;  // source index (1..K) per raw sample
#>   matrix[N_source_raw, J] source_raw;  // raw source isotope measurements
#>   matrix[K, J] tdf_mean_data;  // TDF isotope means, one row per source
#>   matrix[K, J] tdf_sd_data;  // TDF isotope SDs, one row per source
#> }
#> transformed data {
#>   matrix[K, J] tdf_mean = tdf_mean_data;  // alias: TDF means as supplied
#>   matrix[K, J] tdf_sd = tdf_sd_data;  // alias: TDF SDs as supplied
#> }
#> parameters {
#>   simplex[K] p_global;  // global (population-average) source proportions
#>   matrix[2, N_re_Region] z_re_Region;  // std-normal, non-centered
#>   vector<lower=0>[2] sd_re_Region;  // group-level SD, one per (term x ILR dim)
#>   cholesky_factor_corr[2] Lcorr_re_Region;  // group-level correlation (Cholesky factor)
#>   vector<lower=0, upper=20>[J] resid_prop;  // MixSIAR residual-error factor, scales process variance
#>   matrix[K, J] source_mean;  // estimated source isotope means (raw source data)
#>   matrix<lower=0>[K, J] source_sd;  // estimated source isotope SDs (raw source data)
#>   array[K] cholesky_factor_corr[J] Lcorr_source;  // per-source isotope correlation (Cholesky factor)
#> }
#> transformed parameters {
#>   vector[D] ilr_global;  // forward ILR of p_global (Egozcue et al. 2003, eq. 25)
#>   for (d in 1:D) {
#>     ilr_global[d] = sqrt(d / (d + 1.0)) * log(fmax(gmean(p_global, d), 1e-12) / fmax(p_global[d + 1], 1e-12));  // floored to avoid log(0) at a simplex boundary
#>   }
#>   
#>   matrix[N, D] eta = rep_matrix(to_row_vector(ilr_global), N);  // per-mixture-sample ILR predictor: global baseline only (no covariates)
#>   
#>   // group-level term: ((Intercept) | Region)
#>   matrix[2, N_re_Region] b_re_Region = diag_pre_multiply(sd_re_Region, Lcorr_re_Region) * z_re_Region;  // scaled, correlated group-level effects
#>   for (re_i in 1:N) {  // add this term's group-level effect to each mixture sample's eta
#>     for (re_m in 1:M_re_Region) {
#>       for (re_d in 1:D) {
#>         eta[re_i, re_d] += Z_re_Region[re_i, re_m] * b_re_Region[(re_m - 1) * D + re_d, grp_re_Region[re_i]];
#>       }
#>     }
#>   }
#>   
#>   array[N] simplex[K] p;  // source proportions, one simplex per mixture sample
#>   for (i in 1:N) {  // inverse-ILR each mixture sample's eta onto the source simplex
#>     p[i] = softmax(V * to_vector(eta[i]));
#>   }
#>   
#>   matrix[N, J] mu;  // expected mixture isotope value
#>   matrix[N, J] proc_var;  // source/TDF variance propagated into the mixture
#>   for (i in 1:N) {
#>     for (j in 1:J) {  // mu[i, j] = proportion-weighted source + TDF means
#>       real mij = 0;  // expected isotope value
#>       real vij = 1e-8;  // source/TDF variance propagated into the mixture
#>       for (k in 1:K) {  // accumulate each source's contribution, weighted by its proportion
#>         mij += p[i, k] * (source_mean[k, j] + tdf_mean[k, j]);
#>         vij += square(p[i, k]) * (square(source_sd[k, j]) + square(tdf_sd[k, j]));
#>       }
#>       mu[i, j] = mij;
#>       proc_var[i, j] = vij;
#>     }
#>   }
#>   
#>   array[K] matrix[J, J] L_source_cov;  // per-source isotope covariance (Cholesky factor)
#>   array[K] matrix[J, J] source_cov;  // per-source isotope covariance (variance + correlation)
#>   for (k in 1:K) {
#>     L_source_cov[k] = diag_pre_multiply(to_vector(source_sd[k]), Lcorr_source[k]);
#>     source_cov[k] = multiply_lower_tri_self_transpose(L_source_cov[k]);
#>   }
#>   
#>   array[N] matrix[J, J] Omega;  // full source/TDF covariance propagated into the mixture
#>   for (i in 1:N) {
#>     Omega[i] = add_diag(rep_matrix(0, J, J), to_vector(proc_var[i]));  // diagonal: per-isotope process variance
#>     for (j1 in 1:(J - 1)) {
#>       for (j2 in (j1 + 1):J) {
#>         real cij = 0;  // off-diagonal process covariance between isotopes j1, j2
#>         for (k in 1:K) {
#>           cij += square(p[i, k]) * source_cov[k][j1, j2];
#>         }
#>         Omega[i][j1, j2] = cij;
#>         Omega[i][j2, j1] = cij;
#>       }
#>     }
#>   }
#>   
#>   array[N] matrix[J, J] L_Sigma;  // mixture covariance (Cholesky factor)
#>   for (i in 1:N) {
#>     L_Sigma[i] = diag_pre_multiply(sqrt(resid_prop), cholesky_decompose(Omega[i]));  // MixSIAR residual-error factor scales process covariance
#>   }
#> }
#> model {
#>   // prior: global (population-average) source proportions
#>   p_global ~ dirichlet(alpha_dirichlet);
#>   
#>   // priors: group-level effects
#>   to_vector(z_re_Region) ~ std_normal();
#>   sd_re_Region ~ student_t(3, 0, 1);
#>   Lcorr_re_Region ~ lkj_corr_cholesky(1);
#>   
#>   // priors: MixSIAR residual-error factor (scales process variance)
#>   resid_prop[1] ~ uniform(0, 20);  // d13C
#>   resid_prop[2] ~ uniform(0, 20);  // d15N
#>   
#>   // priors + sub-model: raw source data (source- and isotope-specific)
#>   source_mean[1, 1] ~ normal(2.65907, 7.01714);  // Beaver, d13C
#>   source_sd[1, 1] ~ student_t(3, 0, 0.701714);
#>   source_mean[1, 2] ~ normal(-12.3577, 7.73755);  // Beaver, d15N
#>   source_sd[1, 2] ~ student_t(3, 0, 0.773755);
#>   source_mean[2, 1] ~ normal(-7.62958, 7.22181);  // Deer, d13C
#>   source_sd[2, 1] ~ student_t(3, 0, 0.722181);
#>   source_mean[2, 2] ~ normal(12.5269, 6.85126);  // Deer, d15N
#>   source_sd[2, 2] ~ student_t(3, 0, 0.685126);
#>   source_mean[3, 1] ~ normal(12.4633, 7.77402);  // Hare, d13C
#>   source_sd[3, 1] ~ student_t(3, 0, 0.777402);
#>   source_mean[3, 2] ~ normal(1.46555, 3.66305);  // Hare, d15N
#>   source_sd[3, 2] ~ student_t(3, 0, 0.366305);
#>   
#>   // priors: per-source isotope correlation (Cholesky factor)
#>   Lcorr_source[1] ~ lkj_corr_cholesky(1);  // Beaver
#>   Lcorr_source[2] ~ lkj_corr_cholesky(1);  // Deer
#>   Lcorr_source[3] ~ lkj_corr_cholesky(1);  // Hare
#>   
#>   // likelihood for each raw source measurement (isotopes jointly, per-source correlation)
#>   for (n in 1:N_source_raw) {
#>     to_vector(source_raw[n]) ~ multi_normal_cholesky(to_vector(source_mean[source_idx[n]]), L_source_cov[source_idx[n]]);
#>   }
#>   
#>   // mixture likelihood
#>   for (i in 1:N) {  // isotopes jointly, full process (+ residual) covariance
#>     to_vector(y[i]) ~ multi_normal_cholesky(to_vector(mu[i]), L_Sigma[i]);
#>   }
#> }
#> generated quantities {
#>   vector[N] log_lik;  // joint log density per mixture sample (for loo)
#>   matrix[N, J] y_rep;
#>   for (i in 1:N) {
#>     log_lik[i] = multi_normal_cholesky_lpdf(to_vector(y[i]) | to_vector(mu[i]), L_Sigma[i]);
#>     y_rep[i] = to_row_vector(multi_normal_cholesky_rng(to_vector(mu[i]), L_Sigma[i]));
#>   }
#> }
```
