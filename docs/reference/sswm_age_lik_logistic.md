# Logistic continuous-covariate SSWM likelihood

Constructs a negative log-likelihood function for site-level selection
under the strong-selection, weak-mutation (SSWM) model, with selection
intensity modeled as a logistic function of a continuous covariate:

## Usage

``` r
sswm_age_lik_logistic(
  rates_tumors_with,
  rates_tumors_without,
  covariate_data,
  selection_results = NULL
)
```

## Arguments

- rates_tumors_with:

  vector of site-specific mutation rates for all tumors with variant

- rates_tumors_without:

  vector of site-specific mutation rates for all eligible tumors without
  variant

- covariate_data:

  A `data.table` containing `Unique_Patient_Identifier` and
  `covariate_value`. The `covariate_value` column contains the
  sample-level continuous covariate and must be coercible to numeric.

## Value

A function of a three-element numeric parameter vector corresponding to
`L`, `k`, and `m`. The returned function evaluates the negative
log-likelihood for the variant under the logistic continuous-covariate
SSWM model.

## Details

\$\$ \gamma(x) = \frac{\exp(L)}{1 + \exp\[-k(x - m)\]}. \$\$

Here, \\\exp(L)\\ determines the upper asymptote, \\k\\ determines the
steepness or growth/decay rate of the transition, and \\m\\ is the
midpoint of the transition. A positive \\k\\ gives an increasing curve
and a negative \\k\\ gives a decreasing curve.

This function is a likelihood-function factory intended for use with
[`ces_variant_logistic()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_logistic.md).
For each variant, baseline mutation rates and sample information are
supplied automatically by the cancer effect size inference workflow.

## See also

[`ces_variant_logistic`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/ces_variant_logistic.md),
[`sswm_age_lik`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/sswm_age_lik.md),
[`sswm_age_lik_sigmoid`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/sswm_age_lik_sigmoid.md)

Other continuous selection likelihoods:
[`sswm_age_lik()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/sswm_age_lik.md),
[`sswm_age_lik_sigmoid()`](https://townsend-lab-yale.github.io/cancereffectsizeR/reference/sswm_age_lik_sigmoid.md)
